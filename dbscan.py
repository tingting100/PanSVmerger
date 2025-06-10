from collections import defaultdict, Counter
from sklearn.cluster import DBSCAN
import numpy as np
import os
from vcf_parse import read_sv_results, df_to_sv_dict


# k-mer计数函数
def get_kmer_count_vector(seq, k):
    counts = Counter()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        counts[kmer] += 1
    return counts


# 自动计算alpha参数
def auto_alpha_by_length_std(sequences, lambda_param=2.5):
    lengths = [len(seq) for seq in sequences]
    if len(lengths) < 2:
        return 0.5
    std = np.std(lengths)
    mean_len = np.mean(lengths)
    return 1 / (1 + lambda_param * (std / mean_len)) if mean_len != 0 else 0.5


# 计算距离矩阵
def compute_weighted_distance_matrix(sequences, alpha, k):
    print(f"开始计算距离矩阵...")
    n = len(sequences)
    kmer_counts = [get_kmer_count_vector(seq, k) if len(seq) >= k else Counter() for seq in sequences]
    all_kmers = sorted(set(k for c in kmer_counts for k in c))
    kmer_index = {kmer: i for i, kmer in enumerate(all_kmers)}

    freq_matrix = np.zeros((n, len(all_kmers)))
    for i, counts in enumerate(kmer_counts):
        total = sum(counts.values())
        for kmer, cnt in counts.items():
            freq_matrix[i][kmer_index[kmer]] = cnt / total if total != 0 else 0.0

    matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            len_i, len_j = len(sequences[i]), len(sequences[j])
            len_diff = abs(len_i - len_j) / max(len_i, len_j) if max(len_i, len_j) != 0 else 0.0
            i_short, j_short = len_i < k, len_j < k

            if not i_short and not j_short:
                min_sum = np.minimum(freq_matrix[i], freq_matrix[j]).sum()
                max_sum = np.maximum(freq_matrix[i], freq_matrix[j]).sum()
                jaccard = 1.0 - (min_sum / max_sum if max_sum != 0 else 0.0)
            elif i_short and j_short:
                jaccard = 0.0
            else:
                jaccard = 1.0

            matrix[i][j] = matrix[j][i] = alpha * jaccard + (1 - alpha) * len_diff
    return matrix


# 重新分配噪声标签
def reassign_noise_labels(labels):
    print(f"开始重新分配噪声标签...")
    labels = np.array(labels)
    if len(labels) == 0:
        return labels
    max_label = labels[labels != -1].max() if np.any(labels != -1) else -1
    next_label = max_label + 1
    for i in range(len(labels)):
        if labels[i] == -1:
            labels[i] = next_label
            next_label += 1
    return labels


# 单位点处理函数
def process_one_position(chrom, pos, seq_list):
    print(f"开始处理 {chrom} {pos}位点...")
    seq_types = ["REF"] + ["ALT"] * (len(seq_list) - 1)
    filtered_indices = [i for i, seq in enumerate(seq_list) if len(seq) >= 50]
    filtered_seqs = [seq_list[i] for i in filtered_indices]
    result_lines = []

    if not filtered_seqs:
        for i in range(len(seq_list)):
            result_lines.append(f"{chrom}\t{pos}\t{i + 1}\t{seq_types[i]}\t-1\t{seq_list[i]}\tFILTERED")
        return result_lines

    if len(filtered_seqs) == 1:
        for i in range(len(seq_list)):
            if i in filtered_indices:
                status = "SINGLE"
                cluster = 0
            else:
                status = "FILTERED"
                cluster = -1
            result_lines.append(f"{chrom}\t{pos}\t{i + 1}\t{seq_types[i]}\t{cluster}\t{seq_list[i]}\t{status}")
        return result_lines

    alpha = auto_alpha_by_length_std(filtered_seqs)
    dist_matrix = compute_weighted_distance_matrix(filtered_seqs, alpha, k=5)
    dbscan = DBSCAN(eps=0.5, min_samples=2, metric='precomputed')
    labels = dbscan.fit_predict(dist_matrix)
    labels = reassign_noise_labels(labels)

    cluster_map = {filtered_indices[i]: labels[i] for i in range(len(filtered_indices))}

    for i in range(len(seq_list)):
        if i in filtered_indices:
            status = "CLUSTERED"
            cluster = cluster_map[i]
        else:
            status = "FILTERED"
            cluster = -1
        result_lines.append(f"{chrom}\t{pos}\t{i + 1}\t{seq_types[i]}\t{cluster}\t{seq_list[i]}\t{status}")

    return result_lines


# 主聚类函数（串行处理）
def cluster_sv_sequentially(sv_records, output_file):
    print(f"开始串行聚类处理，共{len(sv_records)}个位点...")
    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, 'w') as f_out:
        f_out.write("chrom\tpos\tseq_index\ttype\tcluster\tsequence\tstatus\n")
        processed = 0

        for (chrom, pos), seq_list in sv_records.items():
            if not seq_list:
                continue

            result_lines = process_one_position(chrom, pos, seq_list)
            for line in result_lines:
                f_out.write(line + "\n")

            processed += 1
            if processed % 100 == 0:
                print(f"已处理{processed}个位点...")

    print(f"聚类完成，结果保存至{output_file}")


if __name__ == "__main__":
    vcf_path = "CPC.Phase1.CHM13v2.vcf"
    output_file = "cluster_results_sequential.tsv"

    df = read_sv_results("CPC.Phase1.CHM13v2_SV.tsv")
    if df is None:
        exit(1)

    sv_records = df_to_sv_dict(df)
    cluster_sv_sequentially(sv_records, output_file)