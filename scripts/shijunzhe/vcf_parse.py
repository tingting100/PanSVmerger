from collections import defaultdict
import pandas as pd


def parse_vcf(vcf_path, min_sv_len=50, output_path=None):

    sv_records = defaultdict(list)

    with open(vcf_path, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 5:
                continue

            chrom = fields[0]
            pos = int(fields[1])
            ref = fields[3]
            alt_field = fields[4]

            alt_list = alt_field.split(",")

            # 过滤：只有当ref长度 >= min_sv_len，才加入REF序列
            if len(ref) >= min_sv_len:
                sv_records[(chrom, pos)].append({
                    'type': 'REF',
                    'sequence': ref
                })

            # 添加符合条件的ALT序列，且ALT长度也必须 >= min_sv_len
            for alt_allele in alt_list:
                len_diff = abs(len(ref) - len(alt_allele))
                if len_diff > min_sv_len and len(alt_allele) >= min_sv_len:
                    sv_records[(chrom, pos)].append({
                        'type': 'ALT',
                        'sequence': alt_allele
                    })


    total = sum(len(v) for v in sv_records.values())
    print(f"提取出 {total} 条结构变异（分布在 {len(sv_records)} 个位点），长度差 > {min_sv_len}bp，且REF/ALT序列长度均 >= {min_sv_len}bp。")

    if output_path:
        save_sv_records(sv_records, output_path)

    return sv_records


def save_sv_records(sv_records, output_path):
    with open(output_path, 'w') as f:
        f.write("CHROM\tPOS\tTYPE\tSEQUENCE\n")
        for (chrom, pos) in sorted(sv_records.keys(), key=lambda x: (x[0], x[1])):
            for record in sv_records[(chrom, pos)]:
                f.write(f"{chrom}\t{pos}\t{record['type']}\t{record['sequence']}\n")
    print(f"结果已保存到: {output_path}")


# 读取TSV文件
def read_sv_results(file_path):

    try:
        df = pd.read_csv(file_path, sep='\t')
        print(f"成功读取 {len(df)} 条SV记录")
        return df
    except FileNotFoundError:
        print(f"错误: 文件 {file_path} 不存在")
        return None
    except Exception as e:
        print(f"读取文件时出错: {e}")
        return None


def df_to_sv_dict(df):
    sv_dict = defaultdict(list)
    for _, row in df.iterrows():
        chrom = row['CHROM']
        pos = row['POS']
        seq = row['SEQUENCE']
        sv_dict[(chrom, pos)].append(seq)
    return sv_dict


if __name__ == "__main__":

    vcf_file = "hprc-v1.0-mc-chm13.vcf"
    output_file = "hprc-v1.0-mc-chm13_SV.tsv"

    sv_records = parse_vcf(vcf_file, min_sv_len=50, output_path=output_file)
