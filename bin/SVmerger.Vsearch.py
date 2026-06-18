#!/usr/bin/env python3
"""
PanSVmerger - VSEARCH‑based clustering with length‑based sub‑clustering for multiallelic SVs

Usage:
    python vsearch_clustering.py --input_vcf input.vcf.gz --output_vcf output.vcf [options]

This script runs VSEARCH cluster_fast on all alleles at each locus, then further splits
clusters if sequence lengths differ by >10% (or user-defined threshold).
"""

import os
import sys
import argparse
import gzip
import tempfile
import shutil
import subprocess
from collections import Counter, defaultdict

# --------------------------- Path detection ---------------------------
VSEARCH_PATH = "/public/share/ruanjue/yangting/software/vsearch-2.26.1-linux-x86_64/bin/vsearch"
if not VSEARCH_PATH:
    sys.stderr.write("Warning: 'vsearch' not found in PATH. Set VSEARCH_PATH environment variable.\n")
    VSEARCH_PATH = 'vsearch'  # fallback


# --------------------------- I/O functions ---------------------------
def read_vcf(vcf_file):
    """Read VCF, return variants and header."""
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
    Returns list of allele indices or None if missing.
    """
    if gt_str in ('.', './.', '.|.'):
        return None
    if ploidy == 1:
        if gt_str == '.':
            return None
        try:
            return [int(gt_str)]
        except ValueError:
            return None
    else:
        if '/' in gt_str:
            parts = gt_str.split('/')
        elif '|' in gt_str:
            parts = gt_str.split('|')
        else:
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


# --------------------------- VSEARCH clustering core ---------------------------
def write_locus_fasta(chrom, pos, ref, alts, fasta_path):
    """Write REF and ALT sequences to FASTA file."""
    with open(fasta_path, 'w') as f:
        f.write(f">{chrom}-{pos}-ref\n{ref}\n")
        for j, alt in enumerate(alts, 1):
            f.write(f">{chrom}-{pos}-alt{j}\n{alt}\n")


def run_vsearch(input_fasta, output_uc, threads, identity):
    """Run VSEARCH cluster_fast, return True if successful."""
    cmd = [
        VSEARCH_PATH,
        "--cluster_fast", input_fasta,
        "--id", str(identity),
        "--strand", "both",
        "--uc", output_uc,
        "--threads", str(threads),
        "--maxseqlength", "100000",
        "--minseqlength", "1"   # can be made configurable if needed
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
        return True
    except subprocess.CalledProcessError as e:
        sys.stderr.write(f"VSEARCH error: {e.stderr}\n")
        return False


def parse_vsearch_uc(uc_file):
    """
    Parse UC file, return clusters dict:
        (chrom,pos) -> cluster_id -> list of (alt_id, length, is_ref)
    """
    clusters = defaultdict(lambda: defaultdict(list))
    if not os.path.exists(uc_file) or os.path.getsize(uc_file) == 0:
        return clusters

    with open(uc_file, 'r') as f:
        for line in f:
            if not line.startswith(('S', 'H')):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            try:
                cluster_id = int(fields[1])
                length = int(fields[2])
                seq_id = fields[8]
                # seq_id format: chrom-pos-ref or chrom-pos-altN
                parts = seq_id.split('-')
                if len(parts) < 3:
                    continue
                chrom = parts[0]
                pos = int(parts[1])
                # last part
                last = parts[-1]
                if last == 'ref':
                    alt_id = 0
                    is_ref = True
                elif last.startswith('alt'):
                    alt_id = int(last.replace('alt', ''))
                    is_ref = False
                else:
                    continue
                clusters[(chrom, pos)][cluster_id].append((alt_id, length, is_ref))
            except (ValueError, IndexError):
                continue
    return clusters


def refine_clusters_by_length(clusters, length_cut=0.1):
    """
    For each VSEARCH cluster, further split if length difference > length_cut * previous_length.
    Returns new cluster mapping: (chrom,pos) -> alt_id -> cluster_id
    """
    new_clusters = {}
    for (chrom, pos), cluster_dict in clusters.items():
        # Gather all sequences with their cluster_id
        seq_info = []  # (alt_id, length, is_ref, original_cluster_id)
        for cid, seqs in cluster_dict.items():
            for alt_id, length, is_ref in seqs:
                seq_info.append((alt_id, length, is_ref, cid))

        # Sort by original cluster id, then by length descending
        seq_info.sort(key=lambda x: (x[3], -x[1]))

        # Assign new cluster IDs
        new_cluster_id = 0
        ref_new_cluster = None
        mapping = {}  # alt_id -> new_cluster_id

        # Process each original cluster separately
        current_orig_cid = None
        for alt_id, length, is_ref, orig_cid in seq_info:
            if orig_cid != current_orig_cid:
                # new original cluster
                current_orig_cid = orig_cid
                # start a new new-cluster
                new_cluster_id += 1  # will be adjusted later to make REF=0
                prev_len = None
                # we'll handle within this cluster
                # but we need to process first sequence
                # We'll use a sub-loop
            # We need to re-process by grouping original clusters
        # Better: process each original cluster independently
        final_mapping = {}
        global_new_id = 0
        ref_cluster = None
        # We'll collect all new clusters across original clusters
        for orig_cid, seqs in cluster_dict.items():
            # sort by length descending
            seqs_sorted = sorted(seqs, key=lambda x: x[1], reverse=True)
            subcluster_id = None
            prev_len = None
            for alt_id, length, is_ref in seqs_sorted:
                if subcluster_id is None:
                    subcluster_id = global_new_id
                    global_new_id += 1
                else:
                    # check if length difference exceeds threshold
                    if prev_len - length > length_cut * prev_len:
                        subcluster_id = global_new_id
                        global_new_id += 1
                final_mapping[alt_id] = subcluster_id
                prev_len = length
                if is_ref:
                    ref_cluster = subcluster_id

        # Now remap so REF cluster becomes 0
        if ref_cluster is not None:
            remap = {}
            new_id = 0
            # first assign REF cluster to 0
            remap[ref_cluster] = 0
            # then others
            for cid in sorted(set(final_mapping.values())):
                if cid != ref_cluster:
                    new_id += 1
                    remap[cid] = new_id
            final_mapping = {alt_id: remap[cid] for alt_id, cid in final_mapping.items()}

        new_clusters[(chrom, pos)] = final_mapping

    return new_clusters


def cluster_locus_with_vsearch(fasta_path, threads, identity, length_cut=0.1):
    """
    Run VSEARCH on locus, parse and refine clusters.
    Returns mapping: alt_id -> cluster_id
    """
    # Create temporary UC file
    temp_dir = os.path.dirname(fasta_path)
    uc_file = os.path.join(temp_dir, "vsearch_out.uc")
    if not run_vsearch(fasta_path, uc_file, threads, identity):
        return {}

    clusters = parse_vsearch_uc(uc_file)
    if not clusters:
        return {}

    # Refine by length
    refined = refine_clusters_by_length(clusters, length_cut)
    if (fasta_path is None):  # dummy, but we need to extract chrom,pos
        # but we have it from filename? We'll pass key from caller.
        pass
    # Since refine_clusters_by_length returns dict with (chrom,pos) keys,
    # we need to extract the one for this locus.
    # We'll modify refine_clusters_by_length to return mapping directly for the locus.
    # Better: rewrite refine function to accept a single locus's cluster dict.
    # We'll adjust.
    return refined


# We'll rewrite refine to be per-locus and integrated.
def cluster_locus(fasta_path, threads, identity, length_cut=0.1):
    """Run VSEARCH and return cluster mapping for this locus."""
    temp_dir = os.path.dirname(fasta_path)
    uc_file = os.path.join(temp_dir, "vsearch_out.uc")
    if not run_vsearch(fasta_path, uc_file, threads, identity):
        return {}

    clusters = parse_vsearch_uc(uc_file)
    if not clusters:
        return {}

    # clusters is dict with one key (chrom,pos)
    # We'll extract the first (and only) key
    for (chrom, pos), cluster_dict in clusters.items():
        # refine within this locus
        final_mapping = {}
        global_new_id = 0
        ref_cluster = None
        for orig_cid, seqs in cluster_dict.items():
            seqs_sorted = sorted(seqs, key=lambda x: x[1], reverse=True)
            subcluster_id = None
            prev_len = None
            for alt_id, length, is_ref in seqs_sorted:
                if subcluster_id is None:
                    subcluster_id = global_new_id
                    global_new_id += 1
                else:
                    if prev_len - length > length_cut * prev_len:
                        subcluster_id = global_new_id
                        global_new_id += 1
                final_mapping[alt_id] = subcluster_id
                prev_len = length
                if is_ref:
                    ref_cluster = subcluster_id

        # Remap so REF cluster is 0
        if ref_cluster is not None:
            remap = {ref_cluster: 0}
            new_id = 0
            for cid in sorted(set(final_mapping.values())):
                if cid != ref_cluster:
                    new_id += 1
                    remap[cid] = new_id
            final_mapping = {alt_id: remap[cid] for alt_id, cid in final_mapping.items()}

        return final_mapping

    return {}


# --------------------------- VCF update with ploidy ---------------------------
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
    total_alleles = 0

    for g in genotypes:
        alleles = split_genotype(g, ploidy)
        if alleles is None:
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
                new_alleles.append(str(orig))
                total_alleles += 1

        if ploidy == 1:
            if len(new_alleles) == 1:
                new_genotypes.append(new_alleles[0])
            else:
                new_genotypes.append(','.join(new_alleles))
        else:
            sep = '|' if '|' in g else '/'
            new_genotypes.append(sep.join(new_alleles))

    # Build new ALT list
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
            new_alts.append('N')

    # Update INFO
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
        description="PanSVmerger - VSEARCH-based clustering with length sub-clustering."
    )
    parser.add_argument("--input_vcf", required=True, help="Input VCF (can be .gz)")
    parser.add_argument("--output_vcf", required=True, help="Output VCF after merging")
    parser.add_argument("--new_cluster_file", default="cluster_mapping.txt",
                        help="Output file for cluster mapping (chrom,pos,alt_id,cluster_id)")
    parser.add_argument("--ploidy", type=int, default=2, choices=[1, 2],
                        help="Ploidy: 1 for haploid, 2 for diploid (default: 2)")
    parser.add_argument("--threads", type=int, default=4, help="VSEARCH threads")
    parser.add_argument("--identity", type=float, default=0.9,
                        help="VSEARCH identity threshold for clustering (default: 0.9)")
    parser.add_argument("--length_cut", type=float, default=0.1,
                        help="Proportional length difference to split VSEARCH clusters (default: 0.1)")
    parser.add_argument("--batch_size", type=int, default=1000,
                        help="Progress reporting interval")
    args = parser.parse_args()

    # Read input
    variants, header = read_vcf(args.input_vcf)
    write_vcf_header(args.output_vcf, header)

    total = len(variants)
    all_clusters = {}

    # Create a temporary directory for this run
    temp_dir = tempfile.mkdtemp(prefix="PanSVmerger_vsearch_")
    try:
        for idx, variant in enumerate(variants):
            chrom, pos, vid, ref, alts, info, fmt, genotypes = variant
            if idx % args.batch_size == 0:
                sys.stderr.write(f"Processing {idx}/{total} variants...\n")

            if len(alts) <= 1:
                write_vcf_variant(args.output_vcf, variant)
                all_clusters[(chrom, pos)] = {0: 0}
                if len(alts) == 1:
                    all_clusters[(chrom, pos)][1] = 1
                continue

            # Create temporary FASTA for this locus
            fd, fasta_path = tempfile.mkstemp(suffix='.fasta', dir=temp_dir,
                                              prefix=f"{chrom}_{pos}_")
            os.close(fd)
            write_locus_fasta(chrom, pos, ref, alts, fasta_path)

            try:
                cluster_mapping = cluster_locus(fasta_path, args.threads,
                                                args.identity, args.length_cut)
                if not cluster_mapping:
                    # fallback: treat each ALT as separate
                    cluster_mapping = {0: 0}
                    for i in range(1, len(alts)+1):
                        cluster_mapping[i] = i
                # Ensure all ALT present
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
                write_vcf_variant(args.output_vcf, variant)

            finally:
                if os.path.exists(fasta_path):
                    os.remove(fasta_path)

        # Save cluster mapping
        save_cluster_mapping(all_clusters, args.new_cluster_file)
        sys.stderr.write(f"Done. Output VCF: {args.output_vcf}\n")
        sys.stderr.write(f"Cluster mapping: {args.new_cluster_file}\n")

    finally:
        shutil.rmtree(temp_dir, ignore_errors=True)


if __name__ == "__main__":
    main()
