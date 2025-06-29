import os
import argparse
import logging
from Bio import SeqIO
from concurrent.futures import ThreadPoolExecutor
from tqdm import tqdm
import pandas as pd
import subprocess
from tempfile import NamedTemporaryFile


def parse_fasta_metadata(header):
    parts = header.split("_")
    if len(parts) >= 4:
        chrom = parts[0]  # 如 CHM13v2.chr1
        pos = parts[1]
        type_ = parts[2]
        return chrom, pos, type_
    else:
        return "NA", "NA", "NA"


def run_vsearch(fasta_path, uc_path, identity):
    cmd = [
        "vsearch", "--cluster_fast", fasta_path,
        "--id", str(identity),
        "--uc", uc_path,
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

def process_fasta_file(fasta_path, identity):
    """按50bp bin分组，每组进行vsearch"""
    records = list(SeqIO.parse(fasta_path, "fasta"))
    if not records:
        return []

    chrom, pos, _ = parse_fasta_metadata(records[0].id)

    # 计算序列长度和 length_bin
    data = []
    for r in records:
        length = len(r.seq)
        bin_label = f"{(length//50)*50}-{(length//50)*50+49}"
        data.append((r.id, str(r.seq), bin_label))

    # 分组聚类
    result_rows = []
    df = pd.DataFrame(data, columns=["id", "sequence", "length_bin"])
    for bin_label, bin_df in df.groupby("length_bin"):
        with NamedTemporaryFile(mode="w+", suffix=".fa", delete=False) as tmpf:
            for i, row in bin_df.iterrows():
                tmpf.write(f">{row['id']}\n{row['sequence']}\n")
            tmpf_path = tmpf.name

        uc_file = tmpf_path + ".uc"
        try:
            run_vsearch(tmpf_path, uc_file, identity)
            cluster_map = parse_uc_file(uc_file)
        except Exception as e:
            logging.warning(f"[SKIP] Error running vsearch on bin {bin_label}: {e}")
            cluster_map = {}

        for _, row in bin_df.iterrows():
            cluster = cluster_map.get(row["id"], "NA")
            c, p, t = parse_fasta_metadata(row["id"])
            result_rows.append({
                "chrom": c,
                "pos": p,
                "type": t,
                "length_bin": bin_label,
                "cluster": cluster,
                "sequence": row["sequence"]
            })

        os.remove(tmpf_path)
        if os.path.exists(uc_file):
            os.remove(uc_file)

    return result_rows

def main():
    parser = argparse.ArgumentParser(description="Length-first clustering with vsearch")
    parser.add_argument("--fasta_dir", default="locus_fasta_files", help="Directory with fasta files")
    parser.add_argument("--output", default="length_first_cluster_result.tsv", help="Output TSV file")
    parser.add_argument("--identity", type=float, default=0.9, help="Clustering identity")
    parser.add_argument("--threads", type=int, default=8, help="Number of threads")
    args = parser.parse_args()

    fasta_files = sorted([
        os.path.join(args.fasta_dir, f)
        for f in os.listdir(args.fasta_dir)
        if f.endswith(".fa") and not f.startswith(".")
    ])


    all_results = []
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = [executor.submit(process_fasta_file, f, args.identity) for f in fasta_files]
        for fut in tqdm(futures, desc="Processing files"):
            all_results.extend(fut.result())

    if not all_results:
        logging.warning("❌ 没有成功解析的聚类结果，输出为空。")
        return

    df = pd.DataFrame(all_results)
    df = df[["chrom", "pos", "type", "length_bin", "cluster", "sequence"]]
    df.to_csv(args.output, sep="\t", index=False)
    logging.info(f"✅ 聚类完成，结果已保存到: {args.output}")

if __name__ == "__main__":
    main()
