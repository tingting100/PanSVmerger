#!/bin/bash
# ================================================================
# PanSVmerger - Analysis Pipeline for Pangenome Graph SV Merging
# 
# This script documents the complete analysis workflow for:
#   - Minigraph-Cactus (MC) and PGGB pangenome graphs
#   - Pre-processing, statistics, and three merging strategies
#
# Usage:
#   bash run.sh
#
# Note: Modify the INPUT and PREFIX variables below before running.
# ================================================================

set -e

# ================================================================
# 1. CONFIGURATION (MODIFY THESE)
# ================================================================

# Input data
INPUT_VCF="CHM13v2_HPRCp1_HG002"      # raw graph VCF
PREFIX="MC"                                    # output prefix (MC or PGGB)

# Parameters
PLOIDY=2                                       # 1=haploid, 2=diploid
THREADS=8                                      # threads for parallel processing

# Software paths (adjust to your environment)
PYTHON="~/software/python3.9_env/bin/python3.9"
VSEARCH="~/software/vsearch-2.26.1-linux-x86_64/bin/vsearch"
PERL="perl"
VCFTOOLS="vcftools"
BCFTOOLS="bcftools"

# Reference files
REGION_BED="CHM13v2.telo_cent.complement.bed"

# PanSVmerger scripts
bin="/path/to/PanSVmerger/bin"
GET_SV="$bin/get.SV.pl"
STAT_AC="$bin/stat.countAC.pl"
EXTRACT_AF="$bin/extract_af.py"
SIMILARITY="$bin/allele_allpairs_similarity.py"
JACCARD="$bin/SV.merge.jaccard.py"
VSEARCH_SCRIPT="$bin/SV.merge.Vsearch.py"
LENGTH_SCRIPT="$bin/SV.merge.length.py"

# Chromosomes
CHRS=(chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10
      chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19
      chr20 chr21 chr22 chrX)

# ================================================================
# 2. STEP 1: EXTRACT SV (≥50bp)
# ================================================================
# Filter raw graph VCF to keep only structural variants (≥50bp)

echo "========================================="
echo "[1/6] Extracting SV variants (≥50bp)..."
echo "========================================="

$PERL $GET_SV $INPUT_VCF > ${PREFIX}.raw.SV.vcf

# ================================================================
# 3. STEP 2: AC STATISTICS 
# ================================================================
# Count allele distribution per locus:
#   - How many loci have 1,2,3,...,≥30 ALT alleles?

echo "========================================="
echo "[2/6] Computing AC statistics ..."
echo "========================================="

$PERL $STAT_AC ${PREFIX}.raw.SV.vcf > ${PREFIX}.raw.SV.vcf.AC

# ================================================================
# 4. STEP 3: SPLIT BY CHROMOSOME & FILTER COMPLEX REGIONS
# ================================================================
# Split multi-sample VCF by chromosome and remove:
#   - Telomeres (chr*_telomere*)
#   - Centromeres (chr*_centromere*)
# Using CHM13v2 complement BED file

echo "========================================="
echo "[3/6] Splitting by chromosome and filtering regions..."
echo "========================================="

for chr in "${CHRS[@]}"; do
    echo "  Processing $chr..."
    
    # Extract chromosome
    $VCFTOOLS --vcf ${PREFIX}.raw.SV.vcf --chr $chr \
        --recode --recode-INFO-all --stdout 2>/dev/null \
        > ${PREFIX}.${chr}.vcf
    
    # Remove telomeres and centromeres
    $BCFTOOLS view -T $REGION_BED --threads $THREADS \
        ${PREFIX}.${chr}.vcf \
        -o ${PREFIX}.${chr}.remove_telo_cent.vcf
done

# ================================================================
# 5. STEP 4: EXTRACT AF VALUES 
# ================================================================
# Extract AF values from all ALT alleles across all chromosomes

echo "========================================="
echo "[4/6] Extracting AF values ..."
echo "========================================="

for chr in "${CHRS[@]}"; do
    VCF="${PREFIX}.${chr}.remove_telo_cent.vcf"
    if [ -f "$VCF" ]; then
        echo "  Extracting AF from $chr..."
        $PYTHON $EXTRACT_AF --input_vcf $VCF --output_af ${VCF}.AF
    fi
done

# Combine all AF values
cat ${PREFIX}.*.remove_telo_cent.vcf.AF > ${PREFIX}.all.AF
echo "  Combined AF values saved to ${PREFIX}.all.AF"

# ================================================================
# 6. STEP 5: ALLELE SIMILARITY
# ================================================================
# For each locus, compute all-vs-all sequence similarity using VSEARCH

echo "========================================="
echo "[5/6] Computing allele similarity ..."
echo "========================================="

for chr in "${CHRS[@]}"; do
    VCF="${PREFIX}.${chr}.remove_telo_cent.vcf"
    if [ -f "$VCF" ]; then
        echo "  Computing similarity for $chr..."
        $PYTHON $SIMILARITY \
            --input_vcf $VCF \
            --output_similarity ${PREFIX}.${chr}.similarity.tsv \
            --vsearch_path $VSEARCH \
            --min_alts 2 \
            --threads $THREADS
    fi
done

# ================================================================
# 7. STEP 6: MERGE WITH THREE STRATEGIES
# ================================================================
# Apply three clustering strategies to consolidate redundant alleles:
#   A) Jaccard (k-mer + length) - recommended, best balance
#   B) VSEARCH (global alignment) - highest resolution
#   C) Length (pure length) - fastest, for large datasets

echo "========================================="
echo "[6/6] Merging SVs with three strategies..."
echo "========================================="

for chr in "${CHRS[@]}"; do
    VCF="${PREFIX}.${chr}.remove_telo_cent.vcf"
    [ ! -f "$VCF" ] && continue
    
    echo "  Merging $chr..."
    
    # ---------- Strategy A: Jaccard (k-mer + DBSCAN) ----------
    echo "    A. Jaccard method..."
    $PYTHON $JACCARD \
        --input_vcf $VCF \
        --output_vcf ${PREFIX}.${chr}.Jaccard.vcf \
        --new_cluster_file ${PREFIX}.${chr}.Jaccard.cluster.txt \
        --ploidy $PLOIDY \
        --threads $THREADS \
        --length_threshold 50 \
        --eps_short 0.5 \
        --eps_long 0.3
    
    # ---------- Strategy B: VSEARCH (global alignment) ----------
    echo "    B. VSEARCH method..."
    $PYTHON $VSEARCH_SCRIPT \
        --input_vcf $VCF \
        --output_vcf ${PREFIX}.${chr}.Vsearch.vcf \
        --new_cluster_file ${PREFIX}.${chr}.Vsearch.cluster.txt \
        --ploidy $PLOIDY \
        --threads $THREADS \
        --identity 0.9 \
        --length_cut 0.1
    
    # ---------- Strategy C: Length (pure length) ----------
    echo "    C. Length method..."
    $PYTHON $LENGTH_SCRIPT \
        --input_vcf $VCF \
        --output_vcf ${PREFIX}.${chr}.Length.vcf \
        --new_cluster_file ${PREFIX}.${chr}.Length.cluster.txt \
        --ploidy $PLOIDY \
        --threads $THREADS \
        --eps 0.3
done

# ================================================================
# 8. MERGE CHROMOSOMES & FINAL STATISTICS
# ================================================================
# Concatenate chromosome-level results and compute AC statistics

echo "========================================="
echo "Merging chromosomes and final statistics..."
echo "========================================="

for METHOD in Jaccard Vsearch Length; do
    # Build file list
    VCF_LIST=""
    for chr in "${CHRS[@]}"; do
        VCF="${PREFIX}.${chr}.${METHOD}.vcf"
        [ -f "$VCF" ] && VCF_LIST="$VCF_LIST $VCF"
    done
    
    if [ -n "$VCF_LIST" ]; then
        echo "  Merging $METHOD results..."
        $BCFTOOLS concat -a -o ${PREFIX}.${METHOD}.merged.vcf $VCF_LIST 2>/dev/null
        
        # AC statistics for merged results (Figure 2d)
        $PERL $STAT_AC ${PREFIX}.${METHOD}.merged.vcf > ${PREFIX}.${METHOD}.merged.vcf.AC
    fi
done
