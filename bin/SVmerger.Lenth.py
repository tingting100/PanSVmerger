#!/usr/bin/env python3
"""
PanSVmerger - Length‑based clustering for multiallelic SVs
Usage: python length_clustering.py --input_vcf in.vcf --output_vcf out.vcf [--eps 0.3] [--ploidy 2] [options]
"""

import os
import sys
import argparse
import tempfile
import random
import gzip
from collections import Counter
import numpy as np
from sklearn.cluster import DBSCAN
from Bio import SeqIO


# ---------- I/O functions ----------
def read_vcf(vcf_file):
    """Read VCF, return variants list and header lines."""
    variants = []
    header = []
    if vcf_file.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt'
    else:
        open_func = open
        mode = 'r'

    with open_func(vcf_file, mode) as f:
        for line in f:
            if line.startswith('##'):
                header.append(line.strip().split('\t'))
                continue
            if line.startswith('#CHROM'):
                header.append(line.strip().split('\t'))
                continue
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
        for line_parts in header:
            f.write("\t".join(line_parts) + "\n")


def write_vcf_variant(output_file, variant):
    chrom, pos, vid, ref, alts, info, fmt, genotypes = variant
    with open(output_file, 'a') as f:
        f.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{','.join(alts)}\t.\t.\t{info}\t{fmt}\t"
                + "\t".join(genotypes) + "\n")


def save_cluster_mapping(cluster_mapping, output_file):
    """Save (chrom, pos) -> alt_id -> cluster_id mapping."""
    with open(output_file, 'w') as f:
        for (chrom, pos), alt_dict in sorted(cluster_mapping.items()):
            for alt_id, cluster_id in alt_dict.items():
                f.write(f"{chrom}\t{pos}\t{alt_id}\t{cluster_id}\n")


# ---------- Clustering core ----------
def write_tmp_fasta(chrom, pos, ref, alts, fasta_path):
    """Write all alleles (ref + alts) to a FASTA file."""
    with open(fasta_path, 'w') as f:
        f.write(f">{chrom}-{pos}-ref\n{ref}\n")
        for j, alt in enumerate(alts, 1):
            f.write(f">{chrom}-{pos}-alt{j}\n{alt}\n")


def compute_length_distance_matrix(seqs):
    """Compute pairwise length‑based distance matrix (0..1)."""
    n = len(seqs)
    matrix = np.zeros((n, n))
    for i in range(n):
        li = len(seqs[i])
        for j in range(i + 1, n):
            lj = len(seqs[j])
            if max(li, lj) == 0:
                dist = 0.0
            else:
                dist = abs(li - lj) / max(li, lj)
            matrix[i][j] = matrix[j][i] = dist
    return matrix


def reassign_noise_labels(labels):
    """DBSCAN noise (-1) -> new unique cluster IDs."""
    labels = np.array(labels)
    if len(labels) == 0:
        return labels
    non_noise = labels[labels != -1]
    max_label = non_noise.max() if len(non_noise) > 0 else -1
    next_label = max_label + 1
    for i in range(len(labels)):
        if labels[i] == -1:
            labels[i] = next_label
            next_label += 1
    return labels


def cluster_locus(fasta_path, eps):
    """
    Cluster alleles at one locus using length‑only DBSCAN.
    Returns (cluster_mapping: alt_id->cluster_id, cluster_to_seq: cluster_id->representative_seq)
    """
    try:
        records = list(SeqIO.parse(fasta_path, 'fasta'))
    except Exception:
        print(f"Warning: cannot read FASTA {fasta_path}", file=sys.stderr)
        return {}, {}

    if not records:
        return {}, {}

    ids = [r.id for r in records]
    seqs = [str(r.seq) for r in records]

    # Special cases: 1 or 2 alleles
    if len(seqs) == 1:
        return {0: 0}, {0: seqs[0]}

    if len(seqs) == 2:
        ref_idx = 0
        for i, sid in enumerate(ids):
            if '-ref' in sid:
                ref_idx = i
                break
        alt_idx = 1 - ref_idx
        alt_id = 1
        for sid in ids:
            if '-alt' in sid:
                alt_id = int(sid.split('-')[-1].replace('alt', ''))
                break
        return {0: 0, alt_id: 1}, {0: seqs[ref_idx], 1: seqs[alt_idx]}

    # >= 3 sequences: DBSCAN on length distance matrix
    dist_matrix = compute_length_distance_matrix(seqs)
    db = DBSCAN(eps=eps, min_samples=2, metric='precomputed')
    labels = db.fit_predict(dist_matrix)
    labels = reassign_noise_labels(labels)

    # Find reference index
    ref_idx = 0
    for i, sid in enumerate(ids):
        if '-ref' in sid:
            ref_idx = i
            break

    # Remap labels so that REF cluster becomes 0
    if ref_idx < len(labels):
        ref_label = labels[ref_idx]
        unique_labels = []
        for lab in labels:
            if lab not in unique_labels:
                unique_labels.append(lab)
        if ref_label in unique_labels:
            unique_labels.remove(ref_label)
        ordered = [ref_label] + unique_labels
        remap = {old: new for new, old in enumerate(ordered)}
        labels = np.array([remap[lab] for lab in labels], dtype=int)

    # Build mapping
    cluster_mapping = {}
    cluster_to_seq = {}

    for seq_id, seq, lab in zip(ids, seqs, labels):
        if '-alt' in seq_id:
            alt_id = int(seq_id.split('-')[-1].replace('alt', ''))
            cluster_mapping[alt_id] = lab
            if lab not in cluster_to_seq or len(seq) > len(cluster_to_seq[lab]):
                cluster_to_seq[lab] = seq
        elif '-ref' in seq_id:
            cluster_mapping[0] = lab
            cluster_to_seq[lab] = seq

    return cluster_mapping, cluster_to_seq


# ---------- VCF update with ploidy support ----------
def split_genotype(gt_str, ploidy):
    """
    Split genotype string into list of allele indices based on ploidy.
    Diploid (2): split by '/' or '|', or treat single as 0/0
    Haploid (1): return [int(gt_str)] if not missing
    """
    if gt_str in ('.', './.', '.|.'):
        return None
    if ploidy == 1:
        # Haploid: genotype is a single digit (e.g., "0" or "1")
        if gt_str == '.':
            return None
        try:
            return [int(gt_str)]
        except ValueError:
            return None
    else:
        # Diploid: split by separator
        if '/' in gt_str:
            parts = gt_str.split('/')
        elif '|' in gt_str:
            parts = gt_str.split('|')
        else:
            # Single allele, treat as homozygous
            if gt_str == '.':
                return None
            try:
                return [int(gt_str), int(gt_str)]
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


def update_genotypes_and_ac(variant, all_clusters, cluster_to_seq, ploidy=2):
    chrom, pos, vid, ref, alts, info, fmt, genotypes = variant

    if (chrom, pos) not in all_clusters:
        return variant

    cluster_mapping = all_clusters[(chrom, pos)]

    # Build new genotype mapping: cluster -> new genotype index
    ref_cluster = cluster_mapping[0]
    new_gt_map = {ref_cluster: 0}
    next_gt = 1
    for cluster_id in sorted(set(cluster_mapping.values())):
        if cluster_id != ref_cluster:
            new_gt_map[cluster_id] = next_gt
            next_gt += 1

    # Re‑genotype all samples
    new_genotypes = []
    ac_counts = Counter()
    total_alleles = 0  # AN

    for g in genotypes:
        alleles = split_genotype(g, ploidy)
        if alleles is None:
            new_genotypes.append(g if ploidy == 2 else '.')
            continue

        new_alleles = []
        for orig_gt in alleles:
            if orig_gt is None:
                new_alleles.append('.')
                continue
            if orig_gt in cluster_mapping:
                cid = cluster_mapping[orig_gt]
                new_a = new_gt_map[cid]
                if cid != ref_cluster:
                    ac_counts[new_a] += 1
                new_alleles.append(str(new_a))
                total_alleles += 1
            else:
                # fallback
                new_alleles.append(str(orig_gt))
                total_alleles += 1

        # Format genotype based on ploidy
        if ploidy == 1:
            # Haploid: single value
            if len(new_alleles) == 1:
                new_genotypes.append(new_alleles[0])
            else:
                # Should not happen for haploid
                new_genotypes.append(','.join(new_alleles))
        else:
            # Diploid: join with '/' (unphased) or infer from original
            sep = '/'  # default unphased
            if '|' in g:
                sep = '|'
            new_genotypes.append(sep.join(new_alleles))

    # Build new ALT list from cluster representatives
    new_alts = []
    for cluster_id in sorted(new_gt_map.keys()):
        if cluster_id == ref_cluster:
            continue
        if cluster_id in cluster_to_seq:
            new_alts.append(cluster_to_seq[cluster_id])
        else:
            for alt_id, cid in cluster_mapping.items():
                if cid == cluster_id and alt_id > 0 and alt_id - 1 < len(alts):
                    new_alts.append(alts[alt_id - 1])
                    break

    if not new_alts:
        return variant

    # Update REF sequence
    new_ref = ref
    if ref_cluster in cluster_to_seq:
        new_ref = cluster_to_seq[ref_cluster]

    # Update INFO: AC, AN, AF
    an = total_alleles
    ac_str = ','.join(str(ac_counts[i]) for i in range(1, len(new_alts)+1))
    af_str = ','.join(f"{ac_counts[i]/an:.6f}" if an > 0 else "0" for i in range(1, len(new_alts)+1))

    info_dict = {}
    for field in info.split(';'):
        if '=' in field:
            key, val = field.split('=', 1)
            info_dict[key] = val
        else:
            info_dict[field] = None
    info_dict['AC'] = ac_str
    info_dict['AN'] = str(an)
    info_dict['AF'] = af_str

    new_info = ';'.join(f"{k}={v}" if v is not None else k for k, v in info_dict.items())

    return (chrom, pos, vid, new_ref, new_alts, new_info, fmt, new_genotypes)


# ---------- Main ----------
def main():
    parser = argparse.ArgumentParser(
        description="PanSVmerger - Length‑based multiallelic SV clustering (pure length distance)."
    )
    parser.add_argument("--input_vcf", required=True,
                        help="Input VCF file (can be .gz).")
    parser.add_argument("--output_vcf", required=True,
                        help="Output VCF file after merging.")
    parser.add_argument("--eps", type=float, default=0.3,
                        help="DBSCAN epsilon for length distance (default: 0.3).")
    parser.add_argument("--ploidy", type=int, default=2, choices=[1, 2],
                        help="Ploidy: 1 for haploid, 2 for diploid (default: 2).")
    parser.add_argument("--k", type=int, default=5,
                        help="k-mer size (not used in length-based mode, kept for compatibility).")
    parser.add_argument("--lambda_param", type=float, default=20.0,
                        help="Lambda parameter (not used in length-based mode, kept for compatibility).")
    parser.add_argument("--threads", type=int, default=4,
                        help="Number of threads (not used in length-based mode, kept for compatibility).")
    parser.add_argument("--new_cluster_file", default="cluster_mapping.txt",
                        help="Output file for cluster mapping (chrom,pos,alt_id,cluster_id).")
    args = parser.parse_args()

    variants, header = read_vcf(args.input_vcf)
    write_vcf_header(args.output_vcf, header)

    all_clusters = {}
    all_cluster_seqs = {}

    temp_dir = tempfile.gettempdir()
    pid = os.getpid()
    rand_suffix = random.randint(10000, 99999)
    tmp_fasta = os.path.join(temp_dir, f"tmp_cluster_{pid}_{rand_suffix}.fasta")

    print(f"Processing {len(variants)} variants (ploidy={args.ploidy})...")

    for i, variant in enumerate(variants):
        chrom, pos, vid, ref, alts, info, fmt, genotypes = variant

        if i % 100 == 0:
            print(f"Processed {i}/{len(variants)} variants...")

        if len(alts) == 1:
            write_vcf_variant(args.output_vcf, variant)
            all_clusters[(chrom, pos)] = {0: 0, 1: 1}
            all_cluster_seqs[(chrom, pos)] = {0: ref, 1: alts[0]}
            continue

        write_tmp_fasta(chrom, pos, ref, alts, tmp_fasta)
        cluster_mapping, cluster_to_seq = cluster_locus(tmp_fasta, args.eps)

        if not cluster_mapping:
            write_vcf_variant(args.output_vcf, variant)
            cl_dict = {0: 0}
            seq_dict = {0: ref}
            for alt_idx in range(1, len(alts)+1):
                cl_dict[alt_idx] = alt_idx
                seq_dict[alt_idx] = alts[alt_idx-1]
            all_clusters[(chrom, pos)] = cl_dict
            all_cluster_seqs[(chrom, pos)] = seq_dict
            continue

        max_cluster = max(cluster_mapping.values()) if cluster_mapping else -1
        for alt_idx in range(1, len(alts)+1):
            if alt_idx not in cluster_mapping:
                max_cluster += 1
                cluster_mapping[alt_idx] = max_cluster
                if max_cluster not in cluster_to_seq:
                    cluster_to_seq[max_cluster] = alts[alt_idx-1]

        if 0 not in cluster_mapping:
            cluster_mapping[0] = max_cluster + 1
            cluster_to_seq[cluster_mapping[0]] = ref

        all_clusters[(chrom, pos)] = cluster_mapping
        all_cluster_seqs[(chrom, pos)] = cluster_to_seq

        updated = update_genotypes_and_ac(variant, all_clusters, cluster_to_seq, args.ploidy)
        write_vcf_variant(args.output_vcf, updated)

    if os.path.exists(tmp_fasta):
        os.remove(tmp_fasta)

    save_cluster_mapping(all_clusters, args.new_cluster_file)

    print(f"Done. Output VCF: {args.output_vcf}")
    print(f"Cluster mapping: {args.new_cluster_file}")


if __name__ == "__main__":
    main()
