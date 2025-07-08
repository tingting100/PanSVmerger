import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import os


def group_by_locus_and_create_fasta(tsv_path, output_dir, min_length=50):

    df = pd.read_csv(tsv_path, sep='\t')

    # 标准化列名为大写
    df.columns = df.columns.str.upper()

    os.makedirs(output_dir, exist_ok=True)
    locus_groups = df.groupby(['CHROM', 'POS'])

    total_loci = len(locus_groups)
    valid_loci = 0

    for (chrom, pos), group in locus_groups:
        records = []
        locus_id = f"{chrom}_{pos}"

        for idx, row in group.iterrows():
            seq = str(row['SEQUENCE'])
            seq_type = row['TYPE'].upper()

            if len(seq) >= min_length:
                record_id = f"{locus_id}_{seq_type}_{idx}"
                record = SeqRecord(
                    Seq(seq),
                    id=record_id,
                    description=""
                )
                records.append(record)

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

    locus_fasta_dir = group_by_locus_and_create_fasta(
        "hprc-v1.0-mc-chm13_SV.tsv",
        "hprc_locus_fasta_files",
        min_length=50
    )
