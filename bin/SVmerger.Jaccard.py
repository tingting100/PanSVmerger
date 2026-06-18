#!/usr/bin/env python3
"""
PanSVmerger - Jaccard‑based length‑grouped DBSCAN clustering for multiallelic SVs

Usage:
    python jaccard_clustering.py --input_vcf input.vcf.gz --output_vcf output.vcf [options]

This script clusters alternate alleles at each multiallelic locus using Jaccard distance
and DBSCAN, with separate clustering parameters for short and long sequences.
"""

import os
import sys
import argparse
import gzip
import tempfile
import shutil
from collections import Counter, defaultdict
import numpy as np
from sklearn.cluster import DBSCAN


# --------------------------- I/O functions ---------------------------
def read_vcf(vcf_file):
    """Read VCF file, return list of variants and header lines."""
    variants = []
    header = []
    open_func = gzip.open if vcf_file.endswith('.gz') else open
    mode = 'rt' if vcf_file.endswith('.gz') else 'r'
    with open_func(vcf_file, mode) as f:
        for line in f:
            if line.startswith('##'):
                header.append(line.strip().split('\t'))
            elif line.startswith('#CHROM'):
                header.append(line.strip().split('\t'))
            else:
                fields = line.strip().split('\t')
                chrom, pos, vid = fields[0], int(fields[1]), fields[2]
                ref = fields[3]
                alts = fields[4].split(',') if fields[4] != '.' else []
                info = fields[7]
                fmt = fields[8]
                genotypes = fields[9:]
                variants.append((chrom, pos, vid, ref, alts, info, fmt, genotypes))
    return variants, header


def write_vcf_header(output_file, header):
    with open(output_file, 'w') as f:
        for parts in header:
            f.write("\t".join(parts) + "\n")


def write_vcf_variant(output_file, variant):
    chrom, pos, vid, ref, alts, info, fmt, genotypes = variant
    with open(output_file, 'a') as f:
        f.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{','.join(alts)}\t.\t.\t{info}\t{fmt}\t"
                + "\t".join(genotypes) + "\n")


def save_cluster_mapping(cluster_mapping, output_file):
    """Save (chrom,pos) -> alt_id -> cluster_id mapping."""
    with open(output_file, 'w') as f:
        for (chrom, pos), alt_dict in sorted(cluster_mapping.items()):
            for alt_id, cluster_id in alt_dict.items():
                f.write(f"{chrom}\t{pos}\t{alt_id}\t{cluster_id}\n")


# --------------------------- Genotype helper ---------------------------
def split_genotype(gt_str, ploidy):
    """
    Parse genotype string according to ploidy.
    Returns list of allele indices (integers) or None if missing.
    For diploid: splits by '/' or '|', returns two alleles.
    For haploid: returns a single integer.
    """
    if gt_str in ('.', './.', '.|.'):
        return None
    if ploidy == 1:
        # Haploid: genotype is a single digit
        if gt_str == '.':
            return None
        try:
            return [int(gt_str)]
        except ValueError:
            return None
    else:
        # Diploid
        if '/' in gt_str:
            parts = gt_str.split('/')
        elif '|' in gt_str:
            parts = gt_str.split('|')
        else:
            # Single allele, treat as homozygous
            if gt_str == '.':
                return None
            try:
                val = int(gt_str)
                return [val, val]
            except ValueError:
                return None
        alleles = []
        for p in parts:
            if p == '.':
                alleles.append(None)
            else:
                try:
                    alleles.append(int(p))
                except ValueError:
                    alleles.append(None)
        return alleles


# --------------------------- Jaccard clustering core ---------------------------
def calculate_jaccard_distance(seq1, seq2, k=5):
    """
    Compute Jaccard distance (1 - intersection/union) between two sequences.
    Uses adaptive k-mer size based on sequence length.
    """
    if len(seq1) == 0 or len(seq2) == 0:
        return 1.0
    min_len = min(len(seq1), len(seq2))
    # adaptive k
    if min_len <= 5:
        effective_k = 1
    elif min_len <= 10:
        effective_k = 2
    elif min_len <= 20:
        effective_k = 3
    else:
        effective_k = min(k, 5)
    # Build k-mer sets
    kmers1 = set(seq1[i:i+effective_k] for i in range(len(seq1) - effective_k + 1))
    kmers2 = set(seq2[i:i+effective_k] for i in range(len(seq2) - effective_k + 1))
    if not kmers1 or not kmers2:
        return 1.0
    inter = len(kmers1.intersection(kmers2))
    union = len(kmers1.union(kmers2))
    similarity = inter / union if union > 0 else 0.0
    return 1.0 - similarity


def cluster_by_length_groups(sequences, length_threshold, eps_short, eps_long):
    """
    Cluster sequences using DBSCAN on Jaccard distance, splitting into short/long groups.
    Returns dict {index_in_sequences: cluster_id} with cluster IDs starting from 0.
    """
    n = len(sequences)
    short_indices = [i for i, seq in enumerate(sequences) if len(seq) < length_threshold]
    long_indices = [i for i, seq in enumerate(sequences) if len(seq) >= length_threshold]

    cluster_mapping = {}

    # Cluster short sequences
    if len(short_indices) > 1:
        n_short = len(short_indices)
        dist = np.zeros((n_short, n_short))
        for a in range(n_short):
            for b in range(a+1, n_short):
                d = calculate_jaccard_distance(sequences[short_indices[a]],
                                               sequences[short_indices[b]])
                dist[a, b] = dist[b, a] = d
        db = DBSCAN(eps=eps_short, min_samples=1, metric='precomputed', n_jobs=1)
        labels = db.fit_predict(dist)
        if -1 in labels:
            labels = np.array([0] * n_short)   # treat noise as one cluster
        for idx, orig_idx in enumerate(short_indices):
            cluster_mapping[orig_idx] = labels[idx]
    elif len(short_indices) == 1:
        cluster_mapping[short_indices[0]] = 0

    # Cluster long sequences (offset cluster IDs)
    if len(long_indices) > 1:
        n_long = len(long_indices)
        dist = np.zeros((n_long, n_long))
        for a in range(n_long):
            for b in range(a+1, n_long):
                d = calculate_jaccard_distance(sequences[long_indices[a]],
                                               sequences[long_indices[b]])
                dist[a, b] = dist[b, a] = d
        db = DBSCAN(eps=eps_long, min_samples=1, metric='precomputed', n_jobs=1)
        labels = db.fit_predict(dist)
        if -1 in labels:
            labels = np.array([0] * n_long)
        short_max = max(cluster_mapping.values()) if cluster_mapping else -1
        offset = short_max + 1
        for idx, orig_idx in enumerate(long_indices):
            cluster_mapping[orig_idx] = labels[idx] + offset
    elif len(long_indices) == 1:
        short_max = max(cluster_mapping.values()) if cluster_mapping else -1
        cluster_mapping[long_indices[0]] = short_max + 1

    # Remap cluster IDs to be consecutive
    unique = sorted(set(cluster_mapping.values()))
    remap = {old: new for new, old in enumerate(unique)}
    return {idx: remap[cid] for idx, cid in cluster_mapping.items()}


def cluster_sequences_with_length_groups(fasta_path, length_threshold, eps_short, eps_long):
    """
    Read sequences from FASTA, cluster them, and return mapping from ALT ID to cluster.
    Returns: (id_to_cluster, empty_dict) for compatibility.
    """
    # Read sequences
    sequences = []
    seq_ids = []
    with open(fasta_path, 'r') as f:
        cur_id = None
        cur_seq = []
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if cur_id is not None:
                    sequences.append(''.join(cur_seq))
                    seq_ids.append(cur_id)
                cur_id = line[1:]
                cur_seq = []
            else:
                cur_seq.append(line)
        if cur_id is not None:
            sequences.append(''.join(cur_seq))
            seq_ids.append(cur_id)

    # Find reference
    ref_idx = None
    for i, sid in enumerate(seq_ids):
        if '-ref' in sid:
            ref_idx = i
            break
    if ref_idx is None:
        raise ValueError("No REF sequence found in FASTA")

    cluster_map = cluster_by_length_groups(sequences, length_threshold, eps_short, eps_long)

    # Ensure REF is in cluster 0
    ref_cluster = cluster_map[ref_idx]
    if ref_cluster != 0:
        # Swap: any with 0 becomes ref_cluster, and ref becomes 0
        for idx, cid in cluster_map.items():
            if cid == 0:
                cluster_map[idx] = ref_cluster
        cluster_map[ref_idx] = 0
        # Re-remap to consecutive
        unique = sorted(set(cluster_map.values()))
        remap = {old: new for new, old in enumerate(unique)}
        cluster_map = {idx: remap[cid] for idx, cid in cluster_map.items()}

    # Build alt_id -> cluster mapping
    id_to_cluster = {}
    for idx, sid in enumerate(seq_ids):
        if '-ref' in sid:
            id_to_cluster[0] = cluster_map[idx]
        elif '-alt' in sid:
            # parse alt number
            parts = sid.split('-')
            alt_num = None
            for p in parts:
                if p.startswith('alt'):
                    alt_num = int(p[3:])
                    break
            if alt_num is not None:
                id_to_cluster[alt_num] = cluster_map[idx]

    return id_to_cluster, {}


# --------------------------- VCF update with ploidy support ---------------------------
def ensure_ref_and_alts_in_clusters(variant, clusters):
    """Ensure all ALT alleles have a cluster assignment."""
    chrom, pos, _, _, alts, _, _, _ = variant
    if (chrom, pos) not in clusters:
        clusters[(chrom, pos)] = {}
    if 0 not in clusters[(chrom, pos)]:
        clusters[(chrom, pos)][0] = 0
    max_cid = max(clusters[(chrom, pos)].values(), default=0)
    for i in range(1, len(alts)+1):
        if i not in clusters[(chrom, pos)]:
            max_cid += 1
            clusters[(chrom, pos)][i] = max_cid
    return clusters


def update_genotypes_and_ac(variant, clusters, ploidy=2):
    """
    Update genotypes, AC, AN, AF after clustering, with ploidy support.
    """
    chrom, pos, vid, ref, alts, info, fmt, genotypes = variant
    clusters = ensure_ref_and_alts_in_clusters(variant, clusters)
    cluster_mapping = clusters[(chrom, pos)]

    ref_cluster = cluster_mapping[0]
    cluster_to_gt = {ref_cluster: 0}
    next_gt = 1
    for cid in sorted(set(cluster_mapping.values())):
        if cid != ref_cluster:
            cluster_to_gt[cid] = next_gt
            next_gt += 1

    new_genotypes = []
    ac_counts = Counter()
    total_alleles = 0  # AN

    for g in genotypes:
        alleles = split_genotype(g, ploidy)
        if alleles is None:
            # Missing genotype: keep as is (or '.' for haploid)
            new_genotypes.append(g if ploidy == 2 else '.')
            continue

        new_alleles = []
        for orig in alleles:
            if orig is None:
                new_alleles.append('.')
                continue
            if orig in cluster_mapping:
                cid = cluster_mapping[orig]
                new_a = cluster_to_gt[cid]
                if cid != ref_cluster:
                    ac_counts[new_a] += 1
                new_alleles.append(str(new_a))
                total_alleles += 1
            else:
                # Should not happen
                new_alleles.append(str(orig))
                total_alleles += 1

        # Format genotype according to ploidy
        if ploidy == 1:
            # Haploid: single value
            if len(new_alleles) == 1:
                new_genotypes.append(new_alleles[0])
            else:
                new_genotypes.append(','.join(new_alleles))
        else:
            # Diploid: preserve phasing
            sep = '|' if '|' in g else '/'
            new_genotypes.append(sep.join(new_alleles))

    # Build new ALT list: one per non-ref cluster, using longest representative
    cluster_to_seqs = defaultdict(list)
    for alt_id, cid in cluster_mapping.items():
        if cid != ref_cluster and alt_id > 0 and alt_id - 1 < len(alts):
            cluster_to_seqs[cid].append(alts[alt_id-1])
    new_alts = []
    for cid in sorted(cluster_to_gt.keys()):
        if cid == ref_cluster:
            continue
        if cid in cluster_to_seqs:
            longest = max(cluster_to_seqs[cid], key=len)
            new_alts.append(longest)
        else:
            # Fallback (should not happen)
            new_alts.append('N')

    # Update INFO fields (AC, AN, AF)
    info_dict = {}
    for field in info.split(';'):
        if '=' in field:
            k, v = field.split('=', 1)
            info_dict[k] = v
        else:
            info_dict[field] = None

    an = total_alleles
    ac_str = ','.join(str(ac_counts[i]) for i in range(1, len(new_alts)+1))
    af_str = ','.join(f"{ac_counts[i]/an:.6f}" if an > 0 else "0" for i in range(1, len(new_alts)+1))
    info_dict['AC'] = ac_str
    info_dict['AN'] = str(an)
    info_dict['AF'] = af_str

    new_info = ';'.join(f"{k}={v}" if v is not None else k for k, v in info_dict.items())

    return (chrom, pos, vid, ref, new_alts, new_info, fmt, new_genotypes)


# --------------------------- Main ---------------------------
def main():
    parser = argparse.ArgumentParser(
        description="PanSVmerger - Jaccard‑based length‑grouped DBSCAN clustering."
    )
    parser.add_argument("--input_vcf", required=True,
                        help="Input VCF file (can be .gz)")
    parser.add_argument("--output_vcf", required=True,
                        help="Output VCF file after merging")
    parser.add_argument("--new_cluster_file", default="cluster_mapping.txt",
                        help="Output file for cluster mapping (chrom,pos,alt_id,cluster_id)")
    parser.add_argument("--ploidy", type=int, default=2, choices=[1, 2],
                        help="Ploidy: 1 for haploid, 2 for diploid (default: 2)")
    parser.add_argument("--threads", type=int, default=1,
                        help="Number of threads (not used in DBSCAN, kept for compatibility)")
    parser.add_argument("--length_threshold", type=int, default=50,
                        help="Length threshold to split short/long sequences (bp)")
    parser.add_argument("--eps_short", type=float, default=0.5,
                        help="DBSCAN eps for short sequences (< threshold)")
    parser.add_argument("--eps_long", type=float, default=0.3,
                        help="DBSCAN eps for long sequences (>= threshold)")
    parser.add_argument("--batch_size", type=int, default=1000,
                        help="Progress reporting interval")
    args = parser.parse_args()

    # Check dependencies
    try:
        import numpy as np
        from sklearn.cluster import DBSCAN
    except ImportError as e:
        sys.stderr.write(f"Error: missing dependency: {e}\n")
        sys.stderr.write("Please install: pip install numpy scikit-learn\n")
        sys.exit(1)

    # Read input
    variants, header = read_vcf(args.input_vcf)
    write_vcf_header(args.output_vcf, header)

    total = len(variants)
    all_clusters = {}

    # Create a temporary directory for this run
    temp_dir = tempfile.mkdtemp(prefix="PanSVmerger_jaccard_")
    try:
        for idx, variant in enumerate(variants):
            chrom, pos, vid, ref, alts, info, fmt, genotypes = variant
            if idx % args.batch_size == 0:
                sys.stderr.write(f"Processing {idx}/{total} variants...\n")

            if len(alts) <= 1:
                # No clustering needed
                write_vcf_variant(args.output_vcf, variant)
                all_clusters[(chrom, pos)] = {0: 0}
                if len(alts) == 1:
                    all_clusters[(chrom, pos)][1] = 1
                continue

            # Create a unique temporary FASTA file for this locus
            fd, fasta_path = tempfile.mkstemp(suffix='.fasta', dir=temp_dir,
                                              prefix=f"{chrom}_{pos}_")
            os.close(fd)
            with open(fasta_path, 'w') as f:
                f.write(f">{chrom}-{pos}-ref\n{ref}\n")
                for j, alt in enumerate(alts, 1):
                    f.write(f">{chrom}-{pos}-alt{j}\n{alt}\n")

            try:
                cluster_mapping, _ = cluster_sequences_with_length_groups(
                    fasta_path,
                    args.length_threshold,
                    args.eps_short,
                    args.eps_long
                )

                # Ensure all ALT alleles have a cluster
                for alt_id in range(1, len(alts)+1):
                    if alt_id not in cluster_mapping:
                        max_cid = max(cluster_mapping.values(), default=0)
                        cluster_mapping[alt_id] = max_cid + 1
                if 0 not in cluster_mapping:
                    cluster_mapping[0] = 0

                all_clusters[(chrom, pos)] = cluster_mapping
                updated = update_genotypes_and_ac(variant,
                                                  {(chrom, pos): cluster_mapping},
                                                  args.ploidy)
                write_vcf_variant(args.output_vcf, updated)

            except Exception as e:
                sys.stderr.write(f"Error at {chrom}:{pos}: {e}\n")
                # Write original variant as fallback
                write_vcf_variant(args.output_vcf, variant)

            finally:
                # Remove temporary FASTA file
                if os.path.exists(fasta_path):
                    os.remove(fasta_path)

        # Save cluster mapping
        save_cluster_mapping(all_clusters, args.new_cluster_file)
        sys.stderr.write(f"Done. Output VCF: {args.output_vcf}\n")
        sys.stderr.write(f"Cluster mapping: {args.new_cluster_file}\n")

    finally:
        # Clean up temporary directory
        shutil.rmtree(temp_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
