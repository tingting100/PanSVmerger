#!/usr/bin/env python3.9
"""
SV.merge.jaccard.py
- Parallel DBSCAN clustering
- Short alleles: k-mer similarity only
- Long alleles: Jaccard + length penalty
"""

import os
import argparse
import tempfile
import gzip
import numpy as np
from collections import Counter, defaultdict
from sklearn.cluster import DBSCAN
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm

# =========================
# VCF I/O
# =========================
def read_vcf(path):
    variants, header = [], []
    opener = gzip.open if path.endswith(".gz") else open
    mode = "rt" if path.endswith(".gz") else "r"
    with opener(path, mode) as f:
        for line in f:
            if line.startswith("##"):
                header.append(line.rstrip().split("\t"))
            elif line.startswith("#CHROM"):
                header.append(line.rstrip().split("\t"))
            else:
                fld = line.rstrip().split("\t")
                chrom, pos = fld[0], int(fld[1])
                vid, ref, info, fmt = fld[2], fld[3], fld[7], fld[8]
                alts = [] if fld[4] == "." else fld[4].split(",")
                gts = fld[9:]
                variants.append((chrom, pos, vid, ref, alts, info, fmt, gts))
    return variants, header

def write_tmp_fasta(chrom, pos, ref, alts, path):
    with open(path, "w") as f:
        f.write(f">{chrom}-{pos}-ref\n{ref}\n")
        for i, alt in enumerate(alts, 1):
            f.write(f">{chrom}-{pos}-alt{i}\n{alt}\n")

def write_vcf_header(path, header):
    with open(path, "w") as f:
        for line in header:
            f.write("\t".join(line) + "\n")

def write_vcf_variant(path, var):
    chrom, pos, vid, ref, alts, info, fmt, gts = var
    alts_str = ",".join(alts) if alts else "."
    gts_str = "\t".join(gts)
    with open(path, "a") as f:
        f.write(f"{chrom}\t{pos}\t{vid}\t{ref}\t{alts_str}\t.\t.\t{info}\t{fmt}\t{gts_str}\n")

# =========================
# Cluster mapping (APPEND)
# =========================
def save_cluster_mapping(cluster_mapping, out):
    """Append all loci, never overwrite."""
    with open(out, "a") as f:
        for (chrom, pos), m in cluster_mapping.items():
            for alt_id, cid in m.items():
                f.write(f"{chrom}\t{pos}\t{alt_id}\t{cid}\n")

# =========================
# Distance functions
# =========================
def calculate_distance(seq1, seq2, use_length_penalty=True, alpha=0.7):
    """
    Calculate distance between two sequences.
    - Short: Jaccard only (use_length_penalty=False)
    - Long: Jaccard + length penalty (use_length_penalty=True)
    """
    if not seq1 or not seq2:
        return 1.0

    l1, l2 = len(seq1), len(seq2)
    min_len = min(l1, l2)

    # adaptive k
    if min_len <= 5:
        k = 1
    elif min_len <= 10:
        k = 2
    elif min_len <= 20:
        k = 3
    else:
        k = 5

    s1 = {seq1[i:i+k] for i in range(l1 - k + 1)}
    s2 = {seq2[i:i+k] for i in range(l2 - k + 1)}
    if not s1 or not s2:
        return 1.0

    inter = len(s1 & s2)
    union = len(s1 | s2)
    jaccard = 1.0 - inter / union if union else 1.0

    if not use_length_penalty:
        return jaccard

    length_diff = abs(l1 - l2) / max(l1, l2)
    return alpha * jaccard + (1 - alpha) * length_diff

# =========================
# Clustering core
# =========================
def cluster_by_length_group(sequences, seq_ids, length_threshold, eps_short, eps_long):
    short_idx = [i for i, s in enumerate(sequences) if len(s) < length_threshold]
    long_idx = [i for i, s in enumerate(sequences) if len(s) >= length_threshold]

    cluster = {}

    # ---- Short: k-mer ONLY ----
    if len(short_idx) > 1:
        n = len(short_idx)
        D = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                D[i, j] = D[j, i] = calculate_distance(
                    sequences[short_idx[i]],
                    sequences[short_idx[j]],
                    use_length_penalty=False
                )
        labels = DBSCAN(eps=eps_short, min_samples=1, metric="precomputed").fit_predict(D)
        if -1 in labels:
            labels = np.zeros(n, dtype=int)
        for i, sid in enumerate(short_idx):
            cluster[sid] = labels[i]
    elif len(short_idx) == 1:
        cluster[short_idx[0]] = 0

    # ---- Long: Jaccard + length ----
    if len(long_idx) > 1:
        n = len(long_idx)
        D = np.zeros((n, n))
        for i in range(n):
            for j in range(i + 1, n):
                D[i, j] = D[j, i] = calculate_distance(
                    sequences[long_idx[i]],
                    sequences[long_idx[j]],
                    use_length_penalty=True
                )
        labels = DBSCAN(eps=eps_long, min_samples=1, metric="precomputed").fit_predict(D)
        if -1 in labels:
            labels = np.zeros(n, dtype=int)
        offset = len(set(cluster.values()))
        for i, sid in enumerate(long_idx):
            cluster[sid] = labels[i] + offset
    elif len(long_idx) == 1:
        cluster[long_idx[0]] = len(set(cluster.values()))

    return cluster

def cluster_sequences_with_length_groups(fasta, length_threshold, eps_short, eps_long):
    seqs, ids = [], []
    with open(fasta) as f:
        s, sid = "", ""
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if sid:
                    seqs.append(s)
                    ids.append(sid)
                sid = line[1:]
                s = ""
            else:
                s += line
        if sid:
            seqs.append(s)
            ids.append(sid)

    ref_idx = next(i for i, sid in enumerate(ids) if "ref" in sid.lower())
    cm = cluster_by_length_group(seqs, ids, length_threshold, eps_short, eps_long)

    # ---- Force REF cluster = 0 ----
    ref_cid = cm[ref_idx]
    if ref_cid != 0:
        zero_keys = [k for k, v in cm.items() if v == 0]
        ref_keys = [k for k, v in cm.items() if v == ref_cid]
        for k in zero_keys:
            cm[k] = ref_cid
        for k in ref_keys:
            cm[k] = 0

    uni = sorted(set(cm.values()))
    remap = {o: n for n, o in enumerate(uni)}
    cm = {k: remap[v] for k, v in cm.items()}

    id2c = {}
    for i, sid in enumerate(ids):
        if "ref" in sid.lower():
            id2c[0] = cm[i]
        elif "alt" in sid.lower():
            alt_id = int(sid.split("alt")[-1])
            id2c[alt_id] = cm[i]
    return id2c, {}

# =========================
# VCF update
# =========================
def ensure_ref_and_alts_in_clusters(variant, clusters):
    chrom, pos, _, ref, alts, _, _, _ = variant
    if (chrom, pos) not in clusters:
        clusters[(chrom, pos)] = {}
    if 0 not in clusters[(chrom, pos)]:
        clusters[(chrom, pos)][0] = 0
    mx = max(clusters[(chrom, pos)].values(), default=0)
    for i in range(1, len(alts) + 1):
        if i not in clusters[(chrom, pos)]:
            mx += 1
            clusters[(chrom, pos)][i] = mx
    return clusters

def update_genotypes_and_ac_with_ref_zero(variant, clusters):
    chrom, pos, vid, ref, alts, info, fmt, gts = variant
    clusters = ensure_ref_and_alts_in_clusters(variant, clusters)
    cm = clusters[(chrom, pos)]

    ref_cid = cm[0]
    gt_map = {ref_cid: 0}
    nxt = 1
    for cid in sorted(set(cm.values())):
        if cid != ref_cid:
            gt_map[cid] = nxt
            nxt += 1

    new_gts, ac = [], Counter()
    for gt in gts:
        if gt in (".", "./.", ".|."):
            new_gts.append(gt)
            continue
        sep = "|" if "|" in gt else "/"
        als = gt.split(sep)
        new_a = []
        for a in als:
            if a == ".":
                new_a.append(".")
                continue
            oa = int(a)
            if oa in cm:
                cid = cm[oa]
                ng = gt_map[cid]
                new_a.append(str(ng))
                if ng != 0:
                    ac[ng] += 1
            else:
                new_a.append(a)
        new_gts.append(sep.join(new_a))

    clu_seqs = defaultdict(list)
    for aid, cid in cm.items():
        if aid > 0 and aid - 1 < len(alts):
            clu_seqs[cid].append(alts[aid - 1])

    new_alts = []
    for gid in sorted(gt_map.values()):
        if gid == 0:
            continue
        for cid, g in gt_map.items():
            if g == gid and cid in clu_seqs:
                new_alts.append(max(clu_seqs[cid], key=len))
                break

    # ---- Parse INFO ----
    infod = {}
    for x in info.split(";"):
        if "=" in x:
            k, v = x.split("=", 1)
            infod[k] = v
        else:
            infod[x] = None

    # ---- Calculate AN ----
    an = 0
    for gt in gts:
        if gt not in (".", "./.", ".|."):
            sep = "|" if "|" in gt else "/"
            for a in gt.split(sep):
                if a != ".":
                    an += 1

    # ---- Update AC, AN, AF ----
    infod["AC"] = ",".join(str(ac.get(i, 0)) for i in range(1, len(new_alts) + 1))
    infod["AN"] = str(an)
    af_vals = []
    for i in range(1, len(new_alts) + 1):
        if an > 0:
            af_vals.append(f"{ac.get(i, 0) / an:.6f}")
        else:
            af_vals.append("0.000000")
    infod["AF"] = ",".join(af_vals)

    # ---- Rebuild INFO ----
    new_info = ";".join(f"{k}={v}" if v is not None else k for k, v in infod.items())

    return (chrom, pos, vid, ref, new_alts, new_info, fmt, new_gts)

# =========================
# Main
# =========================
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--input_vcf", required=True)
    ap.add_argument("--output_vcf", required=True)
    ap.add_argument("--new_cluster_file", required=True)
    ap.add_argument("--length_threshold", type=int, default=50)
    ap.add_argument("--eps_short", type=float, default=0.5)
    ap.add_argument("--eps_long", type=float, default=0.3)
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--batch_size", type=int, default=1000)
    args = ap.parse_args()

    # Dependency check
    try:
        import numpy, sklearn
    except ImportError as e:
        raise SystemExit(f"Missing dependency: {e}")

    variants, header = read_vcf(args.input_vcf)
    write_vcf_header(args.output_vcf, header)

    # Init cluster file
    open(args.new_cluster_file, "w").close()

    futures = {}
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        for idx, var in enumerate(variants):
            if idx % args.batch_size == 0:
                print(f"Submitted {idx}/{len(variants)} variants")

            if len(var[4]) <= 1:
                write_vcf_variant(args.output_vcf, var)
                continue

            tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False)
            write_tmp_fasta(var[0], var[1], var[3], var[4], tmp.name)
            tmp.close()

            futures[executor.submit(
                cluster_sequences_with_length_groups,
                tmp.name,
                args.length_threshold,
                args.eps_short,
                args.eps_long
            )] = (var, tmp.name)

        for fut in tqdm(as_completed(futures), total=len(futures), desc="Clustering"):
            var, tmp_path = futures[fut]
            try:
                cm, _ = fut.result()
                save_cluster_mapping({(var[0], var[1]): cm}, args.new_cluster_file)
                upd = update_genotypes_and_ac_with_ref_zero(var, {(var[0], var[1]): cm})
                write_vcf_variant(args.output_vcf, upd)
            except Exception as e:
                print(f"Error at {var[0]}:{var[1]}: {e}")
                write_vcf_variant(args.output_vcf, var)
            finally:
                if os.path.exists(tmp_path):
                    os.remove(tmp_path)

    print("Done.")

if __name__ == "__main__":
    main()

