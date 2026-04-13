import os
import argparse
from collections import Counter, defaultdict
import shutil
import gzip
import re
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
from Bio import SeqIO

def read_vcf(vcf_file):
    # Reads a VCF file and extracts variant information
    variants = []
    header = []
    if vcf_file.endswith('.gz'):
        open_func = gzip.open
        mode = 'rt'
    else:
        open_func = open
        mode = 'r'
    
    with open_func(vcf_file, mode) as file:
        for line in file:
            if line.startswith('##'):
                header.append(line.strip().split('\t'))
                continue
            if line.startswith('#CHROM'):
                header.append(line.strip().split('\t'))
                continue
            fields = line.strip().split('\t')
            chrom, pos, id, ref, info, format_field = fields[0], int(fields[1]), fields[2], fields[3], fields[7], fields[8]
            alts = fields[4].split(',')
            genotypes = fields[9:]
            variants.append((chrom, pos, id, ref, alts, info, format_field, genotypes))
    return variants, header

def write_tmp_fasta(chrom, pos, ref, alts, output_fasta):
    with open(output_fasta, 'w') as file:
        file.write(f">{chrom}-{pos}-ref\n{ref}\n")
        for j, alt in enumerate(alts):
            file.write(f">{chrom}-{pos}-alt{j+1}\n{alt}\n")

def compute_weighted_distance_matrix(sequences, alpha, k):
    """计算距离矩阵，仅基于长度距离"""
    n = len(sequences)
    matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            len_i, len_j = len(sequences[i]), len(sequences[j])
            # 只计算长度距离：归一化的长度差异
            len_diff = abs(len_i - len_j) / max(len_i, len_j) if max(len_i, len_j) != 0 else 0.0
            matrix[i][j] = matrix[j][i] = len_diff
    return matrix

def reassign_noise_labels(labels):
    """重新分配噪声标签（-1）"""
    labels = np.array(labels)
    if len(labels) == 0:
        return labels
    
    # 找出非噪声的最大标签
    non_noise = labels[labels != -1]
    max_label = non_noise.max() if len(non_noise) > 0 else -1
    
    # 为每个噪声分配新标签
    next_label = max_label + 1
    for i in range(len(labels)):
        if labels[i] == -1:
            labels[i] = next_label
            next_label += 1
    
    return labels

def cluster_locus(fasta_path, k, eps, lambda_param):
    """使用基于长度的DBSCAN聚类单个位点"""
    try:
        records = list(SeqIO.parse(fasta_path, 'fasta'))
    except:
        print(f"无法读取FASTA文件: {fasta_path}")
        return {}, {}
    
    if not records:
        return {}, {}
    
    ids = [r.id for r in records]
    seqs = [str(r.seq) for r in records]
    
    # 解析染色体和位置信息
    chrom_pos = os.path.basename(fasta_path).replace('.fa', '').replace('.fasta', '')
    if '-' in chrom_pos:
        chrom, pos = chrom_pos.split('-', 1)
        pos = int(pos)
    else:
        chrom, pos = chrom_pos, 0
    
    # 如果只有一个序列（只有REF）
    if len(seqs) == 1:
        cluster_mapping = {0: 0}  # REF -> cluster 0
        cluster_to_seq = {0: seqs[0]}  # 簇ID -> 序列
        return cluster_mapping, cluster_to_seq
    
    # 如果只有两个序列（REF + 1个ALT）
    if len(seqs) == 2:
        # REF总是簇0
        # ALT是簇1
        cluster_mapping = {}
        cluster_to_seq = {}
        
        # 找到REF和ALT的索引
        ref_idx = 0
        for i, seq_id in enumerate(ids):
            if '-ref' in seq_id:
                ref_idx = i
                break
        
        alt_idx = 1 - ref_idx
        
        # 解析ALT ID
        alt_id = 0
        for seq_id in ids:
            if '-alt' in seq_id:
                alt_part = seq_id.split('-')[-1]
                alt_id = int(alt_part.replace('alt', ''))
                break
        
        cluster_mapping[0] = 0  # REF -> cluster 0
        cluster_mapping[alt_id] = 1  # ALT -> cluster 1
        cluster_to_seq[0] = seqs[ref_idx]
        cluster_to_seq[1] = seqs[alt_idx]
        
        return cluster_mapping, cluster_to_seq
    
    # 3个或更多序列，进行DBSCAN聚类
    # 仅使用长度距离
    dist_matrix = compute_weighted_distance_matrix(seqs, alpha=None, k=None)
    dbscan = DBSCAN(eps=eps, min_samples=2, metric='precomputed')
    labels = dbscan.fit_predict(dist_matrix)
    labels = reassign_noise_labels(labels)
    
    # 找到REF的索引
    ref_idx = 0
    for i, seq_id in enumerate(ids):
        if '-ref' in seq_id:
            ref_idx = i
            break
    
    # 如果找到REF，重映射标签使REF为簇0
    if ref_idx < len(labels):
        ref_label = labels[ref_idx]
        # 收集所有标签
        all_labels = list(dict.fromkeys(labels.tolist()))
        # 确保REF标签在最前面
        if ref_label in all_labels:
            all_labels.remove(ref_label)
        ordered_labels = [ref_label] + all_labels
        remap = {old: new for new, old in enumerate(ordered_labels)}
        labels = np.array([remap[lab] for lab in labels], dtype=int)
    
    # 创建映射关系
    cluster_mapping = {}  # alt_id -> cluster_id
    cluster_to_seq = {}   # cluster_id -> representative_sequence
    
    for seq_id, seq, label in zip(ids, seqs, labels):
        if '-alt' in seq_id:
            # 解析ALT ID
            alt_part = seq_id.split('-')[-1]
            alt_id = int(alt_part.replace('alt', ''))
            cluster_mapping[alt_id] = label
            
            # 记录每个簇的序列（如果同一个簇有多个序列，保留最长的）
            if label not in cluster_to_seq or len(seq) > len(cluster_to_seq[label]):
                cluster_to_seq[label] = seq
        elif '-ref' in seq_id:
            cluster_mapping[0] = label  # REF的alt_id为0
            cluster_to_seq[label] = seq
    
    return cluster_mapping, cluster_to_seq

def ensure_ref_and_alts_in_clusters(variant, clusters):
    chrom, pos, id, ref, alts, info, format_field, genotypes = variant
    
    # Initialize cluster dictionary for the variant's position if it doesn't exist
    if (chrom, pos) not in clusters:
        clusters[(chrom, pos)] = {}
    
    # Ensure the reference sequence (ID = 0) is in the cluster
    if 0 not in clusters[(chrom, pos)]:
        clusters[(chrom, pos)][0] = 0  # Add reference sequence with cluster ID 0

    # Add alternative alleles and ensure they have unique IDs
    for i in range(1, len(alts) + 1):
        if i not in clusters[(chrom, pos)]:
            # Assign a new cluster ID ensuring uniqueness
            clusters[(chrom, pos)][i] = max(clusters[(chrom, pos)].values(), default=0) + 1
    
    return clusters

def update_genotypes_and_ac(variant, clusters, cluster_to_seq):
    chrom, pos, id, ref, alts, info, format_field, genotypes = variant
    
    # 如果没有聚类信息，返回原始variant
    if (chrom, pos) not in clusters:
        return variant
    
    cluster_mapping = clusters[(chrom, pos)]
    
    print(f"Cluster mapping for {chrom}-{pos}: {cluster_mapping}")

    new_genotypes = []
    ac_counts = Counter()

    # Find the reference sequence's cluster ID
    ref_cluster_id = cluster_mapping[0]
    
    # Redefine genotypes based on clusters
    cluster_to_genotype = {ref_cluster_id: 0}
    next_genotype = 1

    for cluster_id in sorted(set(cluster_mapping.values())):
        if cluster_id != ref_cluster_id:
            cluster_to_genotype[cluster_id] = next_genotype
            next_genotype += 1

    for genotype in genotypes:
        if genotype == '.':
            new_genotype = genotype
        else:
            original_genotype = int(genotype)
            if original_genotype in cluster_mapping:
                new_genotype = str(cluster_to_genotype[cluster_mapping[original_genotype]])
                if cluster_mapping[original_genotype] != ref_cluster_id:
                    ac_counts[int(new_genotype)] += 1
            else:
                # Handle singleton sequences
                new_genotype = str(next_genotype)
                next_genotype += 1
                ac_counts[int(new_genotype)] += 1
        new_genotypes.append(new_genotype)

    # Update ALT sequences: take the representative sequence from each cluster
    new_alts = []
    
    # 使用cluster_to_seq中的序列作为代表序列
    for cluster_id in sorted(cluster_to_genotype.keys()):
        if cluster_id != ref_cluster_id:
            if cluster_id in cluster_to_seq:
                new_alts.append(cluster_to_seq[cluster_id])
            else:
                # 如果簇没有代表序列，使用原始ALT中的第一个
                for alt_id, cid in cluster_mapping.items():
                    if cid == cluster_id and alt_id > 0 and alt_id - 1 < len(alts):
                        new_alts.append(alts[alt_id - 1])
                        break
    
    # 如果 new_alts 为空，直接返回原始 variant
    if not new_alts:
        return variant

    # Update INFO field with new AC values
    ac_values = [str(ac_counts[i]) for i in range(1, len(new_alts) + 1)]
    info_fields = info.split(';')
    updated_info_fields = []
    for field in info_fields:
        if field.startswith('AC='):
            updated_info_fields.append('AC=' + ','.join(ac_values))
        else:
            updated_info_fields.append(field)

    # Calculate AF values
    an = sum(1 for genotype in genotypes if genotype != '.')
    # Update INFO field with new AN value
    an_updated = False
    for i, field in enumerate(updated_info_fields):
        if field.startswith('AN='):
            updated_info_fields[i] = f'AN={an}'
            an_updated = True
            break
    if not an_updated:
        updated_info_fields.append(f'AN={an}')

    af_values = [f"{ac_counts[i] / an:.6f}" for i in range(1, len(new_alts) + 1)]
 
    # Update INFO field with new AF values
    af_updated = False
    for i, field in enumerate(updated_info_fields):
        if field.startswith('AF='):
            updated_info_fields[i] = 'AF=' + ','.join(af_values)
            af_updated = True
            break
    if not af_updated:
        updated_info_fields.append('AF=' + ','.join(af_values))
    
    # 更新REF序列
    new_ref = ref
    if 0 in cluster_mapping and cluster_mapping[0] in cluster_to_seq:
        new_ref = cluster_to_seq[cluster_mapping[0]]
    
    updated_info = ';'.join(updated_info_fields)
    
    return (chrom, pos, id, new_ref, new_alts, updated_info, format_field, new_genotypes)

def write_vcf_header(output_file, header):
    with open(output_file, 'w') as file:
        for i in header:
            file.write("\t".join(i) + "\n")

def write_vcf_variant(output_file, variant):
    with open(output_file, 'a') as file:
        chrom, pos, id, ref, alts, info, format_field, genotypes = variant
        file.write("{}\t{}\t{}\t{}\t{}\t.\t.\t{}\t{}\t{}\n".format(
            chrom, pos, id, ref, ",".join(alts), info, format_field, "\t".join(genotypes)))

def save_cluster_mapping(cluster_mapping, output_file):
    # Save the cluster mapping to a file
    with open(output_file, 'w') as file:
        for (chrom, pos), alt_dict in sorted(cluster_mapping.items()):
            for alt_id, cluster_id in alt_dict.items():
                file.write(f"{chrom}\t{pos}\t{alt_id}\t{cluster_id}\n")

def main(input_vcf, output_vcf, tmp_fasta, cluster_prefix, new_cluster_file, threads, eps=0.3, k=5, lambda_param=20.0):
    variants, header = read_vcf(input_vcf)
    write_vcf_header(output_vcf, header)
    
    # 存储所有聚类信息
    all_clusters = {}
    all_cluster_seqs = {}
    
    print(f"Processing {len(variants)} variants...")
    
    for i, variant in enumerate(variants):
        chrom, pos, id, ref, alts, info, format_field, genotypes = variant
        
        if i % 100 == 0:
            print(f"Processed {i}/{len(variants)} variants...")
        
        if len(alts) == 1:
            # 直接输出
            write_vcf_variant(output_vcf, variant)
            # 记录聚类信息（REF为0，ALT为1）
            all_clusters[(chrom, pos)] = {0: 0, 1: 1}
            all_cluster_seqs[(chrom, pos)] = {0: ref, 1: alts[0]}
        else:
            # 创建临时FASTA文件
            write_tmp_fasta(chrom, pos, ref, alts, tmp_fasta)
            
            # 使用DBSCAN聚类
            cluster_mapping, cluster_to_seq = cluster_locus(tmp_fasta, k, eps, lambda_param)
            
            if not cluster_mapping:
                # 如果没有聚类信息，直接输出原始变异
                write_vcf_variant(output_vcf, variant)
                # 为每个ALT分配唯一的簇ID
                cluster_dict = {0: 0}
                seq_dict = {0: ref}
                for alt_idx in range(1, len(alts) + 1):
                    cluster_dict[alt_idx] = alt_idx
                    seq_dict[alt_idx] = alts[alt_idx - 1]
                all_clusters[(chrom, pos)] = cluster_dict
                all_cluster_seqs[(chrom, pos)] = seq_dict
                continue
            
            # 确保所有ALT都在聚类中
            for alt_idx in range(1, len(alts) + 1):
                if alt_idx not in cluster_mapping:
                    # 分配新的簇ID
                    max_cluster = max(cluster_mapping.values()) if cluster_mapping else 0
                    cluster_mapping[alt_idx] = max_cluster + 1
            
            # 确保REF在聚类中
            if 0 not in cluster_mapping:
                cluster_mapping[0] = 0
            
            # 记录聚类信息
            all_clusters[(chrom, pos)] = cluster_mapping
            all_cluster_seqs[(chrom, pos)] = cluster_to_seq
            
            # 更新VCF条目
            updated_variant = update_genotypes_and_ac(variant, all_clusters, cluster_to_seq)
            write_vcf_variant(output_vcf, updated_variant)
    
    # 保存聚类映射
    save_cluster_mapping(all_clusters, new_cluster_file)
    
    print(f"Processing complete! Results saved to {output_vcf}")
    print(f"Cluster mapping saved to {new_cluster_file}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process VCF file and cluster variants using length-based DBSCAN.")
    parser.add_argument("--input_vcf", required=True, help="Input VCF file")
    parser.add_argument("--output_vcf", required=True, help="Output VCF file")
    parser.add_argument("--tmp_fasta", default="tmp_cluster.fasta", help="Temporary FASTA file for clustering")
    parser.add_argument("--cluster_prefix", default="cluster", help="Prefix for cluster output files")
    parser.add_argument("--new_cluster_file", default="cluster_mapping.txt", help="File to save cluster mapping")
    parser.add_argument("--eps", type=float, default=0.3, help="DBSCAN eps parameter")
    parser.add_argument("--k", type=int, default=5, help="k-mer size (not used for length-based clustering)")
    parser.add_argument("--lambda_param", type=float, default=20.0, help="Lambda parameter (not used for length-based clustering)")
    parser.add_argument("--threads", type=int, default=4, help="Threads (kept for compatibility)")
    
    args = parser.parse_args()
    
    main(
        args.input_vcf,
        args.output_vcf,
        args.tmp_fasta,
        args.cluster_prefix,
        args.new_cluster_file,
        args.threads,
        eps=args.eps,
        k=args.k,
        lambda_param=args.lambda_param
    )
