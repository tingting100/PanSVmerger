import os
import glob
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
import numpy as np
from sklearn.cluster import DBSCAN
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm

DBSCAN_EPS         = 0.75
DBSCAN_MIN_SAMPLES = 2
# VAMP_TIMEOUT       = 300

def run_vampire_anno(fa, prefix, threads):
    cmd = ["vampire", "anno", "-t", str(threads), fa, prefix]
    try:
        subprocess.run(cmd,
                       stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL,
                       check=True,
                       # timeout=VAMP_TIMEOUT
                       )
        return True
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        return False

def load_concise(prefix):
    f = prefix + ".concise.tsv"
    if not os.path.exists(f):
        return pd.DataFrame()
    df = pd.read_csv(f, sep="\t", comment="#", low_memory=False)

    if "motif" not in df.columns or not any(c in df.columns for c in ("rep_num","copy_num")):
        return pd.DataFrame()
    col = "rep_num" if "rep_num" in df.columns else "copy_num"
    df["copies"] = pd.to_numeric(df[col], errors="coerce")
    return df.dropna(subset=["motif","copies"])

def cluster_motifs_strict(motifs):

    unique = list(motifs)
    return {m: idx for idx, m in enumerate(unique)}

def dbscan_by_copies(copies):
    X  = copies.reshape(-1,1)
    Xs = StandardScaler().fit_transform(X)
    labs = DBSCAN(eps=DBSCAN_EPS,
                  min_samples=DBSCAN_MIN_SAMPLES)\
           .fit_predict(Xs)
    # 把噪声(-1) 升格为各自独立簇
    if -1 in labs:
        m = labs.max()
        for i,l in enumerate(labs):
            if l == -1:
                m += 1
                labs[i] = m
    return labs

def process_locus(fa_path, threads, out_dir):
    locus     = os.path.basename(fa_path).replace(".fa","")
    locus_dir = os.path.join(out_dir, locus)
    os.makedirs(locus_dir, exist_ok=True)

    prefix = os.path.join(locus_dir, locus)

    # 1) 注释
    if not run_vampire_anno(fa_path, prefix, threads):
        return None

    # 2) 读取 concise.tsv
    df = load_concise(prefix)
    if df.empty:
        return None

    # 3) 第一阶段：strict motif clustering
    motif_map = cluster_motifs_strict(df["motif"].unique())
    df["motif_cluster"] = df["motif"].map(motif_map)

    # 4) 第二阶段：copy-number DBSCAN clustering
    df["copy_cluster"] = -1
    for mc, sub in df.groupby("motif_cluster", sort=False):
        df.loc[sub.index, "copy_cluster"] = dbscan_by_copies(sub["copies"].values)

    # 5) 输出结果
    out_tsv = os.path.join(locus_dir, f"{locus}.clustered.tsv")
    df.to_csv(out_tsv, sep="\t", index=False)
    return locus

def main():
    parser = argparse.ArgumentParser(description="Two-stage TR clustering")
    parser.add_argument("-i","--fasta_dir", required=True,
                        help="Directory of per-locus .fa files")
    parser.add_argument("-t","--threads", type=int, default=1,
                        help="Parallel loci & vampire threads")
    parser.add_argument("-o","--out_dir", default="clustered_results",
                        help="Base output directory")
    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    fastas = sorted(glob.glob(os.path.join(args.fasta_dir, "*.fa")))[:100]

    results = []
    with ThreadPoolExecutor(max_workers=args.threads) as exe:
        futures = {
            exe.submit(process_locus, fa, args.threads, args.out_dir): fa
            for fa in fastas
        }
        for fut in tqdm(as_completed(futures),
                        total=len(futures),
                        desc="Processing loci",
                        unit="locus"):
            try:
                locus = fut.result()
                if locus:
                    results.append(locus)
            except:
                pass

    print(f"\n✅ Done. {len(results)}/{len(fastas)} loci clustered.")
    print("Results in:", args.out_dir)

if __name__ == "__main__":
    main()
