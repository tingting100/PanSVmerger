import os
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
from sklearn.cluster import DBSCAN
from tqdm import tqdm

from vcf_parse import read_sv_results, df_to_sv_dict


def get_kmer_count_vector(seq, k):
    counts = Counter()
    for i in range(len(seq) - k + 1):
        kmer = seq[i:i + k]
        counts[kmer] += 1
    return counts


def auto_alpha_by_length_std(sequences, lambda_param=2.5):
    lengths = [len(seq) for seq in sequences]
    if len(lengths) < 2:
        return 0.5
    std = np.std(lengths)
    mean_len = np.mean(lengths)
    return 1 / (1 + lambda_param * (std / mean_len)) if mean_len != 0 else 0.5


def compute_weighted_distance_matrix(sequences, alpha, k):
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


def reassign_noise_labels(labels):
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


def process_one_position(chrom, pos, seq_list, k, eps, lambda_param):
    seq_types = ["REF"] + ["ALT"] * (len(seq_list) - 1)
    result_lines = []

    if len(seq_list) == 1:
        # 只有一条序列，直接归为一簇0
        result_lines.append(f"{chrom}\t{pos}\t{seq_types[0]}\t0\t{seq_list[0]}")
        return result_lines

    alpha = auto_alpha_by_length_std(seq_list, lambda_param)
    dist_matrix = compute_weighted_distance_matrix(seq_list, alpha, k)
    dbscan = DBSCAN(eps=eps, min_samples=2, metric='precomputed')
    labels = dbscan.fit_predict(dist_matrix)
    labels = reassign_noise_labels(labels)

    for i in range(len(seq_list)):
        cluster = labels[i]
        result_lines.append(f"{chrom}\t{pos}\t{seq_types[i]}\t{cluster}\t{seq_list[i]}")
    return result_lines


def cluster_sv_parallel(sv_records, output_file, k, eps, lambda_param, max_workers=8):

    result_lines = []

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {
            executor.submit(process_one_position, chrom, pos, seq_list, k, eps, lambda_param): (chrom, pos)
            for (chrom, pos), seq_list in sv_records.items() if seq_list
        }

        for future in tqdm(as_completed(futures), total=len(futures), desc="DBSCAN clustering"):
            result = future.result()
            result_lines.extend(result)

    with open(output_file, 'w') as f_out:
        f_out.write("chrom\tpos\ttype\tcluster\tsequence\n")
        for line in result_lines:
            f_out.write(line + "\n")


    print(f"\n✓ DBSCAN 聚类完成，结果保存至：{output_file}")


if __name__ == "__main__":
    # 硬编码参数配置区
    tsv_path = "CPC.Phase1.CHM13v2_SV.tsv"
    output_file = "weighted_joint_cluster_result.tsv"
    k = 5
    eps = 0.5
    lambda_param = 20
    max_workers = 4

    df = read_sv_results(tsv_path)
    if df is None:
        exit(1)

    sv_records = df_to_sv_dict(df)

    cluster_sv_parallel(
        sv_records=sv_records,
        output_file=output_file,
        k=k,
        eps=eps,
        lambda_param=lambda_param,
        max_workers=max_workers,
    )
