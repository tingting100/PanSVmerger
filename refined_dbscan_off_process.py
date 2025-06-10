import matplotlib
matplotlib.use('Agg')  # ✅ 修复 RuntimeError: main thread is not in main loop
import pandas as pd
import matplotlib.pyplot as plt
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
from collections import Counter

# ----------- 可配置参数 -----------
tsv_path = 'cluster_results_sequential.tsv'
output_stats_file = 'cluster_similarity_stats.tsv'
output_plot_dir = 'cluster_stats_parallel_threads'
num_threads = 6
k = 5  # k-mer 长度
# ----------------------------------

os.makedirs(output_plot_dir, exist_ok=True)


def get_kmers(seq, k):
    return set(seq[i:i + k] for i in range(len(seq) - k + 1)) if len(seq) >= k else set()

def jaccard(seq1, seq2):
    kmer1 = get_kmers(seq1, k)
    kmer2 = get_kmers(seq2, k)
    if not kmer1 or not kmer2:
        return 0.0
    return len(kmer1 & kmer2) / len(kmer1 | kmer2)

def process_cluster(group_key, df):
    chrom, pos, cluster_id = group_key
    print(f"正在处理 {chrom} {pos}......")
    records = []

    # 找最长序列为中心
    sequences = df['sequence'].tolist()
    center_idx = df['sequence'].str.len().idxmax()
    center_seq = df.loc[center_idx, 'sequence']

    similarities = [jaccard(center_seq, seq) for seq in sequences]
    lengths = [len(seq) for seq in sequences]

    for idx, (i, row) in enumerate(df.iterrows()):
        print(f"正在写入数据......")
        records.append({
            'chrom': row['chrom'],
            'pos': row['pos'],
            'seq_index': row['seq_index'],
            'type': row['type'],
            'cluster': row['cluster'],
            'length': len(row['sequence']),
            'similarity': similarities[idx]  # ✅ 使用 enumerate 保证同步
        })

    # 画图
    print(f"正在绘制图像......")
    length_counter = Counter(lengths)
    plt.figure(figsize=(6, 4))
    plt.bar(length_counter.keys(), length_counter.values())
    plt.xlabel('Sequence Length')
    plt.ylabel('Count')
    plt.title(f'Length Distribution: {chrom}:{pos}, Cluster {cluster_id}')
    plt.tight_layout()
    plt.savefig(os.path.join(output_plot_dir, f'{chrom}_{pos}_cluster{cluster_id}.png'))
    plt.close()

    return records

def main():
    print("🔍 正在读取聚类结果...")
    df = pd.read_csv(tsv_path, sep='\t')

    print("🔎 正在过滤 cluster == -1 的序列...")
    df = df[df['cluster'] != -1].copy()
    df['key'] = list(zip(df['chrom'], df['pos'], df['cluster']))
    groups = list(df.groupby('key'))

    all_records = []
    print(f"🚀 使用 {num_threads} 个线程处理 {len(groups)} 个簇...")

    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(process_cluster, key, group) for key, group in groups]
        for f in tqdm(futures, desc="Processing Clusters"):
            all_records.extend(f.result())

    # 写入 TSV
    print(f"💾 写入聚类统计结果到 {output_stats_file}")
    pd.DataFrame(all_records).to_csv(output_stats_file, sep='\t', index=False)

    print(f"\n✅ 所有簇已处理完毕，结果保存在: {output_stats_file} 和 {output_plot_dir}/")


if __name__ == '__main__':
    main()