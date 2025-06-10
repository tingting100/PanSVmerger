import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os
from glob import glob


def group_by_locus_and_create_fasta(tsv_path, output_dir, min_length=50):
    """
    按位点分组序列并为每个位点创建独立的FASTA文件

    参数:
        tsv_path: 输入TSV文件路径
        output_dir: 输出FASTA文件目录
        min_length: 最小序列长度阈值
    """
    df = pd.read_csv(tsv_path, sep='\t')
    os.makedirs(output_dir, exist_ok=True)
    locus_groups = df.groupby(['CHROM', 'POS'])

    total_loci = len(locus_groups)
    valid_loci = 0

    for (chrom, pos), group in locus_groups:
        records = []
        locus_id = f"{chrom}_{pos}"

        for idx, row in group.iterrows():
            ref = row['REF']
            alt = row['ALT']

            # 过滤短序列
            if len(ref) >= min_length:
                ref_id = f"{locus_id}_REF_{idx}"
                records.append(SeqRecord(
                    Seq(ref),
                    id=ref_id,
                    description=f"locus={locus_id};type=REF;length={len(ref)}"
                ))

            if len(alt) >= min_length:
                alt_id = f"{locus_id}_ALT_{idx}"
                records.append(SeqRecord(
                    Seq(alt),
                    id=alt_id,
                    description=f"locus={locus_id};type=ALT;length={len(alt)}"
                ))

        # 仅保存有有效序列的位点
        if records:
            fasta_path = os.path.join(output_dir, f"{locus_id}.fa")
            SeqIO.write(records, fasta_path, "fasta")
            valid_loci += 1

    print(f"按位点分组完成：")
    print(f"  总位点: {total_loci}")
    print(f"  有效位点(≥1条序列): {valid_loci}")
    print(f"  输出目录: {output_dir}")
    return output_dir


if __name__ == '__main__':
    # 使用示例
    locus_fasta_dir = group_by_locus_and_create_fasta(
        "CPC.Phase1.CHM13v2_SV.tsv",
        "locus_fasta_files",
        min_length=50
    )