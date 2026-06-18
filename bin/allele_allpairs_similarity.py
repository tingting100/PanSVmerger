cat allele_allpairs_similarity.py 
#!/usr/bin/env python3
"""
Compute pairwise sequence similarity for each multiallelic VCF locus using VSEARCH --allpairs_global.

Usage:
    python vcf_allpairs_similarity.py --input_vcf input.vcf.gz --output_similarity all_pairs.tsv
    [--vsearch_path /path/to/vsearch] [--threads 4] [--min_alts 2] [--sample_loci 1000]

For each locus with >= min_alts ALT alleles, all sequences (REF + ALTs) are compared all-vs-all.
Results are written as a TSV file with columns: CHROM  POS  QUERY  TARGET  IDENTITY  ALN_LEN  MISMATCHES  GAPS.
"""

import os
import sys
import argparse
import gzip
import tempfile
import subprocess
import shutil
import random

def read_vcf(vcf_file):
    """Yield (chrom, pos, ref, alts) for each variant, skipping header."""
    open_func = gzip.open if vcf_file.endswith('.gz') else open
    mode = 'rt' if vcf_file.endswith('.gz') else 'r'
    with open_func(vcf_file, mode) as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue
            chrom = fields[0]
            pos = fields[1]
            ref = fields[3]
            alts = fields[4].split(',') if fields[4] != '.' else []
            yield chrom, pos, ref, alts

def write_locus_fasta(chrom, pos, ref, alts, fasta_path):
    """Write sequences for a single locus to FASTA file."""
    with open(fasta_path, 'w') as f:
        f.write(f">{chrom}-{pos}-ref\n{ref}\n")
        for i, alt in enumerate(alts, 1):
            f.write(f">{chrom}-{pos}-alt{i}\n{alt}\n")

def run_vsearch_allpairs(fasta_path, vsearch_path, threads):
    """
    Run vsearch --allpairs_global on fasta_file.
    Returns a list of tuples: (query, target, identity, aln_len, mismatches, gaps)
    """
    cmd = [
        vsearch_path,
        "--allpairs_global", fasta_path,
        "--acceptall",
        "--self",
        "--threads", str(threads),
        "--blast6out", "STDOUT"
    ]
    try:
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        out, err = proc.communicate()
        if proc.returncode != 0:
            sys.stderr.write(f"VSEARCH error for {fasta_path}: {err}\n")
            return []
    except FileNotFoundError:
        sys.stderr.write(f"VSEARCH not found at {vsearch_path}. Please check path.\n")
        return []

    results = []
    for line in out.splitlines():
        parts = line.strip().split('\t')
        if len(parts) < 6:
            continue
        query, target = parts[0], parts[1]
        identity = float(parts[2])   # percent identity
        aln_len = int(parts[3])
        mismatches = int(parts[4])
        gaps = int(parts[5])
        results.append((query, target, identity, aln_len, mismatches, gaps))
    return results

def parse_seq_id(seq_id):
    """Parse chrom, pos, and allele type from sequence ID like 'chr1-12345-ref' or 'chr1-12345-alt1'."""
    parts = seq_id.split('-')
    if len(parts) < 3:
        return None, None, None
    chrom = parts[0]
    pos = parts[1]
    allele = parts[2]  # 'ref' or 'alt1', 'alt2', ...
    return chrom, pos, allele

def main():
    parser = argparse.ArgumentParser(
        description="Compute all-vs-all sequence similarity per VCF locus using VSEARCH."
    )
    parser.add_argument("--input_vcf", required=True, help="Input VCF file (.vcf or .vcf.gz)")
    parser.add_argument("--output_similarity", required=True, help="Output TSV file for pairwise similarities")
    parser.add_argument("--vsearch_path", default="vsearch", help="Path to VSEARCH executable (default: vsearch)")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads (default: 4)")
    parser.add_argument("--min_alts", type=int, default=2, help="Minimum number of ALT alleles to process (default: 2)")
    parser.add_argument("--sample_loci", type=int, default=0,
                        help="If >0, randomly sample this many loci (useful for large VCFs)")
    parser.add_argument("--include_ref", action="store_true",
                        help="Include REF sequence in all-vs-all (default: only ALTs? Actually we always include REF, but this flag toggles? We'll include REF always; but we can add option to exclude REF if needed, but for similarity density we usually want REF too? We'll keep REF always.")
    args = parser.parse_args()

    # Check VSEARCH availability
    if not shutil.which(args.vsearch_path):
        sys.stderr.write(f"Error: VSEARCH not found at '{args.vsearch_path}'. Please install or specify path.\n")
        sys.exit(1)

    # Read all loci (or sample)
    loci = []
    for chrom, pos, ref, alts in read_vcf(args.input_vcf):
        if len(alts) >= args.min_alts:
            loci.append((chrom, pos, ref, alts))
    
    total_loci = len(loci)
    if args.sample_loci > 0 and args.sample_loci < total_loci:
        random.seed(42)
        loci = random.sample(loci, args.sample_loci)
        sys.stderr.write(f"Sampled {args.sample_loci} loci from {total_loci}\n")
    else:
        sys.stderr.write(f"Processing {total_loci} loci with >= {args.min_alts} ALTs\n")

    # Prepare output file
    out_fh = open(args.output_similarity, 'w')
    out_fh.write("CHROM\tPOS\tQUERY\tTARGET\tIDENTITY\tALN_LEN\tMISMATCHES\tGAPS\n")

    # Create temporary directory
    temp_dir = tempfile.mkdtemp(prefix="vcf_similarity_")
    try:
        for idx, (chrom, pos, ref, alts) in enumerate(loci):
            if idx % 100 == 0:
                sys.stderr.write(f"Processed {idx}/{len(loci)} loci...\n")
            
            # Write FASTA
            fd, fasta_path = tempfile.mkstemp(suffix='.fasta', dir=temp_dir)
            os.close(fd)
            write_locus_fasta(chrom, pos, ref, alts, fasta_path)

            # Run VSEARCH
            results = run_vsearch_allpairs(fasta_path, args.vsearch_path, args.threads)
            
            # Write results for this locus
            for query, target, identity, aln_len, mismatches, gaps in results:
                # Extract chrom,pos from query (they should be same)
                # We can also parse from query ID to ensure consistency
                out_fh.write(f"{chrom}\t{pos}\t{query}\t{target}\t{identity:.2f}\t{aln_len}\t{mismatches}\t{gaps}\n")
            
            # Cleanup
            os.remove(fasta_path)

    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)
        out_fh.close()

    sys.stderr.write(f"Done. Results written to {args.output_similarity}\n")

if __name__ == "__main__":
    main()
