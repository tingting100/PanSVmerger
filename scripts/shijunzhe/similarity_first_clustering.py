import os
import argparse
import pandas as pd
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
import subprocess
from collections import defaultdict
from datasketch import MinHash, MinHashLSH
from sourmash import MinHash as SourmashMinHash
from sklearn.cluster import AgglomerativeClustering
import numpy as np
import screed
import re

# MiniHash 配置
K = 5
NUM_PERM = 128
MINIHASH_THRESHOLD = 0.9

# Sourmash 参数
SOURMASH_K = 5
SOURMASH_SCALED = 1000
SOURMASH_THRESHOLD = 0.9

# CD-HIT 配置
CDHIT_EXEC = "cd-hit-est"
CDHIT_IDENTITY = 0.9
CDHIT_OUT_DIR = "cdhit_cluster_tmp"
os.makedirs(CDHIT_OUT_DIR, exist_ok=True)

def run_vsearch(fasta_file, output_uc, identity=0.9):
    cmd = [
        "vsearch", "--cluster_fast", fasta_file,
        "--id", str(identity),
        "--uc", output_uc,
        "--threads", "1", "--quiet"
    ]
    subprocess.run(cmd, check=True)

def parse_uc_file(uc_file):
    cluster_map = {}
    for line in open(uc_file):
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if parts[0] in ("S", "H"):
            cluster_id = parts[1]
            name = parts[8]
            cluster_map[name] = cluster_id
    return cluster_map

def initial_vsearch_clustering(fasta_file, identity=0.9):
    chrom, pos = os.path.basename(fasta_file).replace(".fa", "").split("_", 1)
    uc_file = fasta_file + ".uc"
    try:
        run_vsearch(fasta_file, uc_file, identity)
        cluster_map = parse_uc_file(uc_file)
        records = list(SeqIO.parse(fasta_file, "fasta"))
        return [
            (chrom.replace("CHM13v2.", ""), pos, r.id.split("_")[2], cluster_map.get(r.id, "NA"), str(r.seq))
            for r in records
        ]
    except Exception as e:
        print(f"[ERROR] VSEARCH failed on {fasta_file}: {e}")
        return []

def create_minhash(seq):
    mh = MinHash(num_perm=NUM_PERM)
    for i in range(len(seq) - K + 1):
        mh.update(seq[i:i+K].encode('utf-8'))
    return mh

def initial_minhash_clustering(fasta_file):
    chrom, pos = os.path.basename(fasta_file).replace(".fa", "").split("_", 1)
    records = list(SeqIO.parse(fasta_file, "fasta"))
    if not records:
        return []

    minhashes = {}
    seq_info = {}
    for r in records:
        minhashes[r.id] = create_minhash(str(r.seq))
        seq_info[r.id] = str(r.seq)

    lsh = MinHashLSH(threshold=MINIHASH_THRESHOLD, num_perm=NUM_PERM)
    for sid, mh in minhashes.items():
        lsh.insert(sid, mh)

    visited = set()
    clusters = []
    for sid in minhashes:
        if sid in visited:
            continue
        bucket = list(set(lsh.query(minhashes[sid])) - visited)
        if bucket:
            visited.update(bucket)
            clusters.append(bucket)

    result = []
    for i, cluster in enumerate(clusters):
        for sid in cluster:
            type_ = sid.split("_")[2]
            result.append((chrom.replace("CHM13v2.", ""), pos, type_, i, seq_info[sid]))
    return result

def initial_sourmash_clustering(fasta_file):
    chrom, pos = os.path.basename(fasta_file).replace(".fa", "").split("_", 1)
    seqs, ids, mhs = [], [], []
    for record in screed.open(fasta_file):
        seq = record.sequence.upper()
        mh = SourmashMinHash(n=0, ksize=SOURMASH_K, scaled=SOURMASH_SCALED)
        mh.add_sequence(seq, force=True)
        seqs.append(seq)
        ids.append(record.name)
        mhs.append(mh)

    if len(mhs) <= 1:
        type_ = ids[0].split("_")[2] if mhs else "NA"
        return [(chrom.replace("CHM13v2.", ""), pos, type_, 0, seqs[0])] if mhs else []

    dist_mat = np.zeros((len(mhs), len(mhs)))
    for i in range(len(mhs)):
        for j in range(i+1, len(mhs)):
            dist = 1 - mhs[i].jaccard(mhs[j])
            dist_mat[i][j] = dist_mat[j][i] = dist

    clustering = AgglomerativeClustering(
        metric='precomputed', linkage='average',
        distance_threshold=1 - SOURMASH_THRESHOLD, n_clusters=None
    )
    labels = clustering.fit_predict(dist_mat)

    result = []
    for i, lbl in enumerate(labels):
        type_ = ids[i].split("_")[2]
        result.append((chrom.replace("CHM13v2.", ""), pos, type_, int(lbl), seqs[i]))
    return result

def run_cdhit(fasta_file, identity=0.9):
    basename = os.path.basename(fasta_file).replace(".fa", "")
    output_prefix = os.path.join(CDHIT_OUT_DIR, basename)
    clstr_file = output_prefix + ".clstr"
    cmd = [
        CDHIT_EXEC,
        "-i", fasta_file,
        "-o", output_prefix,
        "-c", str(identity),
        "-n", "5",
        "-T", "1",
        "-M", "8000",
        "-d", "0"
    ]
    subprocess.run(cmd, check=True)
    return clstr_file

def parse_clstr(clstr_file, fasta_file):
    seq_dict = {r.name: str(r.sequence).upper() for r in SeqIO.parse(fasta_file, "fasta")}
    records = []
    cluster_id = -1
    with open(clstr_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">Cluster"):
                cluster_id += 1
            else:
                m = re.search(r">(\S+)\.\.\.", line)
                if m:
                    full_id = m.group(1)
                    if full_id in seq_dict:
                        chrom_part, pos = full_id.split("_", 1)
                        sv_type = "ALT" if "ALT" in full_id else "REF"
                        records.append((chrom_part.replace("CHM13v2.", ""), pos, sv_type, cluster_id, seq_dict[full_id]))
    return records

def initial_cdhit_clustering(fasta_file, identity=0.9):
    try:
        clstr_file = run_cdhit(fasta_file, identity)
        return parse_clstr(clstr_file, fasta_file)
    except Exception as e:
        print(f"[ERROR] CD-HIT failed on {fasta_file}: {e}")
        return []

def recluster_by_length(input_tsv, output_tsv, identity=0.9, threads=8):
    df = pd.read_csv(input_tsv, sep="\t")
    df["length"] = df["sequence"].str.len()
    df["length_bin"] = df["length"].apply(lambda l: f"{(l//50)*50}-{(l//50)*50+49}")

    grouped = df.groupby(["chrom", "pos", "cluster", "length_bin"])
    results = []

    def cluster_group(key_df):
        from tempfile import NamedTemporaryFile
        with NamedTemporaryFile(mode='w+', suffix='.fa', delete=False) as tmp:
            for i, row in key_df.iterrows():
                tmp.write(f">{i}\n{row['sequence']}\n")
            tmp_path = tmp.name

        uc_file = tmp_path + ".uc"
        try:
            run_vsearch(tmp_path, uc_file, identity)
            cluster_map = parse_uc_file(uc_file)
        except:
            cluster_map = {}

        final_rows = []
        for i, row in key_df.iterrows():
            cid = cluster_map.get(str(i), "NA")
            final_rows.append({
                "chrom": row["chrom"],
                "pos": row["pos"],
                "type": row["type"],
                "cluster": row["cluster"],
                "length_bin": row["length_bin"],
                "subcluster": cid,
                "sequence": row["sequence"]
            })
        return final_rows

    with ThreadPoolExecutor(max_workers=threads) as executor:
        futures = [executor.submit(cluster_group, group.copy()) for _, group in grouped]
        for fut in tqdm(futures, desc="Reclustering by length"):
            results.extend(fut.result())

    out_df = pd.DataFrame(results)
    out_df.to_csv(output_tsv, sep="\t", index=False)
    print(f"✓ Reclustering done. Output saved to: {output_tsv}")

def main():
    parser = argparse.ArgumentParser(description="similarity_first SV clustering")
    parser.add_argument("--fasta_dir", default="locus_fasta_files", help="Directory containing FASTA files")
    parser.add_argument("--identity", type=float, default=0.9, help="Sequence identity for clustering")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads")
    parser.add_argument("--output", default="similarity_first_cluster_result.tsv", help="Output TSV file")
    parser.add_argument("--step", choices=["initial", "recluster"], default="initial", help="Step type")
    parser.add_argument("--method", choices=["vsearch", "minihash", "sourmash", "cdhit"], default="vsearch", help="Clustering method")
    args = parser.parse_args()

    if args.step == "initial":
        fasta_files = sorted([
            os.path.join(args.fasta_dir, f)
            for f in os.listdir(args.fasta_dir) if f.endswith(".fa") and not f.startswith(".")
        ])
        all_data = []

        with ThreadPoolExecutor(max_workers=args.threads) as executor:
            if args.method == "vsearch":
                futures = [executor.submit(initial_vsearch_clustering, f, args.identity) for f in fasta_files]
            elif args.method == "minihash":
                futures = [executor.submit(initial_minhash_clustering, f) for f in fasta_files]
            elif args.method == "sourmash":
                futures = [executor.submit(initial_sourmash_clustering, f) for f in fasta_files]
            elif args.method == "cdhit":
                futures = [executor.submit(initial_cdhit_clustering, f, args.identity) for f in fasta_files]

            for fut in tqdm(futures, desc=f"Initial clustering ({args.method})"):
                all_data.extend(fut.result())

        df = pd.DataFrame(all_data, columns=["chrom", "pos", "type", "cluster", "sequence"])
        df.to_csv(args.output, sep="\t", index=False)
        print(f"✓ Initial clustering done. Output saved to: {args.output}")

    elif args.step == "recluster":
        recluster_by_length(args.output, args.output.replace(".tsv", "_reclustered.tsv"),
                            identity=args.identity, threads=args.threads)

if __name__ == "__main__":
    main()
