import os
import argparse
import numpy as np
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
from tqdm import tqdm
from sklearn.cluster import DBSCAN
import pysam

def read_vcf(vcf_file):
    """Read VCF file and return variants list and header"""
    variants = []
    vcf = pysam.VariantFile(vcf_file)
    
    header = vcf.header
    
    for record in vcf:
        variants.append((
            record.chrom,
            record.pos,
            record.ref,
            record.alts if record.alts else [],
            record
        ))
    
    vcf.close()
    return variants, header

def write_tmp_fasta(chrom, pos, ref, alts, tmp_filepath):
    """Write sequences to temporary FASTA file"""
    os.makedirs(os.path.dirname(tmp_filepath), exist_ok=True)
    with open(tmp_filepath, 'w') as f:
        # Write REF sequence
        f.write(f">REF_{chrom}_{pos}\n")
        f.write(f"{ref}\n")
        
        # Write ALT sequences
        for i, alt in enumerate(alts):
            f.write(f">ALT{i}_{chrom}_{pos}\n")
            f.write(f"{alt}\n")

def get_kmer_count_vector(seq, k):
    """Get kmer count vector for a sequence"""
    counts = Counter()
    for i in range(len(seq) - k + 1):
        counts[seq[i:i + k]] += 1
    return counts

def auto_alpha_by_length_std(sequences, lambda_param=2.5):
    """Automatically calculate alpha parameter based on length standard deviation"""
    lengths = [len(seq) for seq in sequences]
    if len(lengths) < 2:
        return 0.5
    std = np.std(lengths)
    mean_len = np.mean(lengths)
    return 1 / (1 + lambda_param * (std / mean_len)) if mean_len != 0 else 0.5

def compute_weighted_distance_matrix(sequences, alpha, k):
    """Compute weighted distance matrix"""
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
                long_seq = sequences[j] if i_short else sequences[i]
                jaccard = 0.0 if len(long_seq) < 50 else 1.0

            matrix[i][j] = matrix[j][i] = alpha * jaccard + (1 - alpha) * len_diff
    return matrix

def reassign_noise_labels(labels):
    """Reassign noise labels"""
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

def process_one_position(chrom, pos, ref, alts, k, eps, lambda_param, tmp_filepath):
    """Process one genomic position using temporary FASTA file"""
    # 写入序列到临时文件
    write_tmp_fasta(chrom, pos, ref, alts, tmp_filepath)
    
    # 从临时文件读取序列
    sequences = []
    seq_ids = []
    
    with open(tmp_filepath) as f:
        current_seq = ''
        for line in f:
            if line.startswith('>'):
                if current_seq:
                    sequences.append(current_seq)
                seq_ids.append(line.strip()[1:])
                current_seq = ''
            else:
                current_seq += line.strip()
        if current_seq:
            sequences.append(current_seq)
    
    # 执行聚类
    alpha = auto_alpha_by_length_std(sequences, lambda_param)
    dist_matrix = compute_weighted_distance_matrix(sequences, alpha, k)
    
    # 处理特殊情况：所有序列太短或距离矩阵全零
    if np.all(dist_matrix == 0):
        # 为每个序列分配唯一标签
        labels = list(range(len(sequences)))
    else:
        dbscan = DBSCAN(eps=eps, min_samples=1, metric='precomputed')
        labels = dbscan.fit_predict(dist_matrix)
        labels = reassign_noise_labels(labels)
    
    # 确保REF序列在cluster 0
    ref_id = f"REF_{chrom}_{pos}"
    ref_idx = seq_ids.index(ref_id)
    ref_label = labels[ref_idx]
    if ref_label != 0:
        # 交换标签
        labels[labels == 0] = ref_label
        labels[labels == ref_label] = 0
    
    # 构建聚类结果
    clusters = {}
    for seq_id, label in zip(seq_ids, labels):
        clusters[seq_id] = int(label)
    
    return {
        "chrom": chrom,
        "pos": pos,
        "clusters": clusters
    }

def update_genotypes(variants, clustering_results):
    """Update genotype information based on clustering"""
    updated_variants = []
    
    for chrom, pos, ref, alts, record in variants:
        if (chrom, pos) not in clustering_results:
            updated_variants.append(record)
            continue
        
        clusters = clustering_results[(chrom, pos)]["clusters"]
        ref_cluster = clusters.get(f"REF_{chrom}_{pos}", 0)
        
        unique_clusters = set(clusters.values())
        if len(unique_clusters) <= 1:
            # 只有一种等位基因，跳过该位点
            updated_variants.append(record)
            continue
 
        # Build cluster to genotype mapping
        cluster_to_gt = {ref_cluster: 0}
        next_gt = 1
        
        # Process ALTs
        new_alts = []
        alt_to_gt = {}
        cluster_representatives = {}

        for i, alt in enumerate(alts):
            alt_key = f"ALT{i}_{chrom}_{pos}"
            cluster_id = clusters.get(alt_key, ref_cluster)
            
            if cluster_id not in cluster_to_gt:
                if cluster_id != ref_cluster:  # New cluster gets new genotype
                    cluster_to_gt[cluster_id] = next_gt
                    next_gt += 1
            
            # Record ALT to genotype mapping
            gt = cluster_to_gt[cluster_id]
            alt_to_gt[i] = gt
            
            # Keep only representative ALTs (first ALT per cluster)
            if gt != 0 and cluster_id not in cluster_representatives:
                cluster_representatives[cluster_id] = alt
                new_alts.append(alt)
        
        # Update genotypes
        for sample in record.samples:
            gt = record.samples[sample]["GT"]
            if gt is None:
                continue
            
            new_gt = []
            for allele in gt:
                if allele is None:
                    new_gt.append(None)
                elif allele == 0:
                    new_gt.append(0)
                else:
                    new_gt.append(alt_to_gt.get(allele - 1, 0))
            
            record.samples[sample]["GT"] = tuple(new_gt)
        
        # Update ALT list
        record.alts = tuple(new_alts)
        updated_variants.append(record)
    
    return updated_variants

def cleanup_temp_dirs(tmp_dir):
    """Clean up empty chromosome directories after all processing is done"""
    if os.path.exists(tmp_dir):
        # 遍历临时目录中的所有染色体目录
        for chrom_dir in os.listdir(tmp_dir):
            chrom_dir_path = os.path.join(tmp_dir, chrom_dir)
            if os.path.isdir(chrom_dir_path):
                # 检查染色体目录是否为空
                if not os.listdir(chrom_dir_path):
                    try:
                        os.rmdir(chrom_dir_path)
                        print(f"Removed empty directory: {chrom_dir_path}")
                    except OSError as e:
                        print(f"Warning: Could not remove directory {chrom_dir_path}: {e}")
                else:
                    print(f"Directory {chrom_dir_path} is not empty, keeping it.")
        
        # 检查主临时目录是否为空
        if not os.listdir(tmp_dir):
            try:
                os.rmdir(tmp_dir)
                print(f"Removed empty temporary directory: {tmp_dir}")
            except OSError as e:
                print(f"Warning: Could not remove temporary directory {tmp_dir}: {e}")
        else:
            print(f"Temporary directory {tmp_dir} is not empty, keeping it.")

def write_output_vcf(output_file, header, variants):
    """Write the output VCF file"""
    with pysam.VariantFile(output_file, 'w', header=header) as out_vcf:
        for variant in variants:
            out_vcf.write(variant)


def main():
    parser = argparse.ArgumentParser(description="SV sequence clustering and genotype updating")
    parser.add_argument("--input", required=True, help="Input VCF file")
    parser.add_argument("--output", required=True, help="Output VCF file")
    parser.add_argument("--k", type=int, default=5, help="k-mer size")
    parser.add_argument("--eps", type=float, default=0.5, help="DBSCAN eps parameter")
    parser.add_argument("--lambda", type=float, default=2.5, dest="lambda_param", help="Lambda parameter for alpha calculation")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    parser.add_argument("--tmp_dir", default="tmp_cluster_files", help="Temporary dir for clustering (optional)")

    args = parser.parse_args()

    # 创建临时目录
    os.makedirs(args.tmp_dir, exist_ok=True)

    # Read VCF
    variants, header = read_vcf(args.input)
    
    # Pre-filter variants
    variants_to_process = []
    clustering_results = {}
    
    # First pass: handle simple cases (single ALT) directly
    for variant in variants:
        chrom, pos, ref, alts, record = variant
        if len(alts) <= 1:
            clusters = {"REF": 0}
            if len(alts) == 1:
                clusters[f"ALT0_{chrom}_{pos}"] = 1
            clustering_results[(chrom, pos)] = {
                "chrom": chrom,
                "pos": pos,
                "clusters": clusters
            }
        else:
            variants_to_process.append((chrom, pos, ref, alts, record))

    # Second pass: process variants with multiple ALTs
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = {}
        for chrom, pos, ref, alts, record in variants_to_process:
            # 生成临时文件名：tmp_dir/chr_pos.fa
            tmp_filename = os.path.join(args.tmp_dir, f"{chrom}_{pos}.fa")
            
            future = executor.submit(
                process_one_position, 
                chrom, pos, ref, alts, 
                args.k, args.eps, args.lambda_param,
                tmp_filename  # 传递临时文件路径
            )
            futures[future] = (chrom, pos, tmp_filename)
        
        # 处理结果并清理临时文件
        for future in tqdm(as_completed(futures), total=len(futures), desc="Processing variants"):
            chrom, pos, tmp_filename = futures[future]
            try:
                result = future.result()
                clustering_results[(chrom, pos)] = result
            except Exception as e:
                print(f"Error processing {chrom}:{pos}: {e}")
                # 出错时返回默认聚类结果
                clusters = {"REF": 0}
                for i in range(len(alts)):
                    clusters[f"ALT{i}_{chrom}_{pos}"] = i + 1
                clustering_results[(chrom, pos)] = {
                    "chrom": chrom,
                    "pos": pos,
                    "clusters": clusters
                }
            finally:
                # 删除临时文件
                if os.path.exists(tmp_filename):
                    try:
                        os.remove(tmp_filename)
                    except Exception as e:
                        print(f"Warning: Could not delete {tmp_filename}: {e}")

    # 清理临时目录（如果为空）
    try:
        if not os.listdir(args.tmp_dir):
            os.rmdir(args.tmp_dir)
            print(f"Temporary directory '{args.tmp_dir}' removed.")
    except OSError:
        print(f"Temporary directory '{args.tmp_dir}' is not empty, keeping it.")
    
    # Update genotypes and write output
    print("Updating genotypes...")
    updated_variants = update_genotypes(variants, clustering_results)
    
    print("Writing output VCF...")
    write_output_vcf(args.output, header, updated_variants)
    print("Processing completed!")

if __name__ == "__main__":
    main()

