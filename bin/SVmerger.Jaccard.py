import os
import argparse
from collections import Counter, defaultdict
import gzip
import numpy as np
from sklearn.cluster import DBSCAN

def read_vcf(vcf_file):
    """读取VCF文件"""
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
    """写入临时FASTA文件"""
    with open(output_fasta, 'w') as file:
        file.write(f">{chrom}-{pos}-ref\n{ref}\n")
        for j, alt in enumerate(alts):
            file.write(f">{chrom}-{pos}-alt{j+1}\n{alt}\n")

def calculate_distance(seq1, seq2, k=5):
    """计算两个序列的Jaccard距离（0-1之间）"""
    if len(seq1) == 0 or len(seq2) == 0:
        return 1.0
    
    len1, len2 = len(seq1), len(seq2)
    min_len = min(len1, len2)
    
    # 自适应k值
    if min_len <= 5:
        effective_k = 1  # 超短序列用 1-mer
    elif min_len <= 10:
        effective_k = 2  # 短序列用 2-mer
    elif min_len <= 20:
        effective_k = 3  # 中等序列用 3-mer
    else:
        effective_k = min(k, 5)  # 长序列用 5-mer 
    # 计算k-mer集合
    kmers1 = set(seq1[i:i+effective_k] for i in range(len1 - effective_k + 1))
    kmers2 = set(seq2[i:i+effective_k] for i in range(len2 - effective_k + 1))
    
    if not kmers1 or not kmers2:
        return 1.0
    
    intersection = len(kmers1.intersection(kmers2))
    union = len(kmers1.union(kmers2))
    
    similarity = intersection / union if union > 0 else 0.0
    return 1.0 - similarity

def cluster_by_length_group(sequences, seq_ids, length_threshold=50, eps_short=0.5, eps_long=0.3):
    """
    根据序列长度分组聚类
    Returns: id_to_cluster mapping {seq_index: cluster_id}
    """
    # 按长度分组
    short_indices = [i for i, seq in enumerate(sequences) if len(seq) < length_threshold]
    long_indices = [i for i, seq in enumerate(sequences) if len(seq) >= length_threshold]
    
    print(f"DEBUG: Short sequences (<{length_threshold}bp): {len(short_indices)}")
    print(f"DEBUG: Long sequences (≥{length_threshold}bp): {len(long_indices)}")
    
    # 初始化cluster mapping
    cluster_mapping = {}
    
    # 聚类短序列
    if len(short_indices) > 1:
        print(f"DEBUG: Clustering short sequences with eps={eps_short}")
        # 构建距离矩阵
        n_short = len(short_indices)
        dist_matrix_short = np.zeros((n_short, n_short))
        
        for idx_i, i in enumerate(short_indices):
            for idx_j, j in enumerate(short_indices):
                if idx_i < idx_j:
                    dist = calculate_distance(sequences[i], sequences[j])
                    dist_matrix_short[idx_i][idx_j] = dist
                    dist_matrix_short[idx_j][idx_i] = dist
        
        # DBSCAN聚类
        dbscan_short = DBSCAN(eps=eps_short, min_samples=1, metric='precomputed', n_jobs=1)
        labels_short = dbscan_short.fit_predict(dist_matrix_short)
        
        print(f"DEBUG: Short sequence labels: {labels_short}")
        
        # 处理噪声点（全部归为一个簇）
        if -1 in labels_short:
            labels_short = np.array([0] * n_short)
        
        # 保存短序列聚类结果
        for idx, seq_idx in enumerate(short_indices):
            cluster_mapping[seq_idx] = labels_short[idx]
    elif len(short_indices) == 1:
        cluster_mapping[short_indices[0]] = 0
    
    # 聚类长序列
    if len(long_indices) > 1:
        print(f"DEBUG: Clustering long sequences with eps={eps_long}")
        # 构建距离矩阵
        n_long = len(long_indices)
        dist_matrix_long = np.zeros((n_long, n_long))
        
        for idx_i, i in enumerate(long_indices):
            for idx_j, j in enumerate(long_indices):
                if idx_i < idx_j:
                    dist = calculate_distance(sequences[i], sequences[j])
                    dist_matrix_long[idx_i][idx_j] = dist
                    dist_matrix_long[idx_j][idx_i] = dist
        
        # DBSCAN聚类
        dbscan_long = DBSCAN(eps=eps_long, min_samples=1, metric='precomputed', n_jobs=1)
        labels_long = dbscan_long.fit_predict(dist_matrix_long)
        
        print(f"DEBUG: Long sequence labels: {labels_long}")
        
        # 处理噪声点（全部归为一个簇）
        if -1 in labels_long:
            labels_long = np.array([0] * n_long)
        
        # 保存长序列聚类结果（偏移cluster id，避免与短序列冲突）
        short_cluster_count = len(set(cluster_mapping.values())) if cluster_mapping else 0
        for idx, seq_idx in enumerate(long_indices):
            cluster_mapping[seq_idx] = labels_long[idx] + short_cluster_count
    elif len(long_indices) == 1:
        # 确定长序列的起始cluster id
        if cluster_mapping:
            next_cluster_id = max(cluster_mapping.values()) + 1
        else:
            next_cluster_id = 0
        cluster_mapping[long_indices[0]] = next_cluster_id
    
    return cluster_mapping

def cluster_sequences_with_length_groups(fasta_path, length_threshold=50, eps_short=0.5, eps_long=0.3):
    """
    根据长度分组聚类序列
    """
    # 读取序列
    sequences = []
    seq_ids = []
    
    with open(fasta_path, 'r') as f:
        current_seq = ""
        current_id = ""
        
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences.append(current_seq)
                    seq_ids.append(current_id)
                current_id = line[1:]
                current_seq = ""
            else:
                current_seq += line
        
        if current_id:
            sequences.append(current_seq)
            seq_ids.append(current_id)
    
    print(f"DEBUG: Total sequences: {len(sequences)}")
    
    # 检查REF序列
    ref_index = -1
    for i, seq_id in enumerate(seq_ids):
        if 'ref' in seq_id.lower():
            ref_index = i
            break
    
    if ref_index == -1:
        raise ValueError("No REF sequence found")
    
    print(f"DEBUG: REF length: {len(sequences[ref_index])}bp")
    
    # 根据长度分组聚类
    cluster_mapping = cluster_by_length_group(sequences, seq_ids, length_threshold, eps_short, eps_long)
    
    # 确保REF所在的簇为簇0
    ref_cluster_id = cluster_mapping[ref_index]
    
    if ref_cluster_id != 0:
        # 交换簇ID：将REF的簇设为0，原来的0簇设为ref_cluster_id
        old_zero_indices = [idx for idx, cid in cluster_mapping.items() if cid == 0]
        ref_indices = [idx for idx, cid in cluster_mapping.items() if cid == ref_cluster_id]
        
        # 交换
        for idx in old_zero_indices:
            cluster_mapping[idx] = ref_cluster_id
        for idx in ref_indices:
            cluster_mapping[idx] = 0
    
    # 重新映射簇ID使其连续
    unique_clusters = sorted(set(cluster_mapping.values()))
    cluster_remap = {old: new for new, old in enumerate(unique_clusters)}
    
    final_mapping = {}
    for idx, cluster_id in cluster_mapping.items():
        final_mapping[idx] = cluster_remap[cluster_id]
    
    # 创建ID到簇的映射
    id_to_cluster = {}
    for idx, seq_id in enumerate(seq_ids):
        if 'ref' in seq_id.lower():
            id_to_cluster[0] = final_mapping[idx]
        elif 'alt' in seq_id.lower():
            parts = seq_id.split('-')
            for part in parts:
                if part.startswith('alt'):
                    try:
                        alt_id = int(part[3:])
                        id_to_cluster[alt_id] = final_mapping[idx]
                    except:
                        alt_id = len([k for k in id_to_cluster.keys() if k > 0]) + 1
                        id_to_cluster[alt_id] = final_mapping[idx]
                    break
    
    # 分析聚类结果
    print("DEBUG: Cluster summary:")
    for cluster_id in sorted(set(final_mapping.values())):
        cluster_indices = [i for i, cid in final_mapping.items() if cid == cluster_id]
        seq_lengths = [len(sequences[i]) for i in cluster_indices]
        seq_types = []
        for i in cluster_indices:
            if i == ref_index:
                seq_types.append("REF")
            elif len(sequences[i]) < length_threshold:
                seq_types.append("short")
            else:
                seq_types.append("long")
        
        print(f"  Cluster {cluster_id}: {len(cluster_indices)} seqs, lengths {seq_lengths}, types {seq_types}")
    
    return id_to_cluster, {}

def ensure_ref_and_alts_in_clusters(variant, clusters):
    chrom, pos, id, ref, alts, info, format_field, genotypes = variant
    
    if (chrom, pos) not in clusters:
        clusters[(chrom, pos)] = {}
    
    # 确保REF在聚类中，且簇ID为0
    if 0 not in clusters[(chrom, pos)]:
        clusters[(chrom, pos)][0] = 0

    # 确保所有ALT在聚类中
    for i in range(1, len(alts) + 1):
        if i not in clusters[(chrom, pos)]:
            clusters[(chrom, pos)][i] = max(clusters[(chrom, pos)].values(), default=0) + 1
    
    return clusters

def update_genotypes_and_ac_with_ref_zero(variant, clusters):
    """
    更新基因型和AC，确保REF的基因型为0
    """
    chrom, pos, id, ref, alts, info, format_field, genotypes = variant
    clusters = ensure_ref_and_alts_in_clusters(variant, clusters)
    cluster_mapping = clusters[(chrom, pos)]

    print(f"=== Processing {chrom}:{pos} ===")
    print(f"REF: {len(ref)}bp")
    print(f"ALTs: {[f'{len(alt)}bp' for alt in alts]}")
    print(f"Cluster mapping: {cluster_mapping}")

    new_genotypes = []
    ac_counts = Counter()

    # === 关键：确保REF的簇映射到基因型0 ===
    # REF的簇ID
    ref_cluster_id = cluster_mapping[0]
    
    # 创建簇ID到基因型的映射
    # REF的簇总是映射到基因型0
    cluster_to_genotype = {ref_cluster_id: 0}
    
    # 其他簇从1开始编号
    next_genotype = 1
    for cluster_id in sorted(set(cluster_mapping.values())):
        if cluster_id != ref_cluster_id:
            cluster_to_genotype[cluster_id] = next_genotype
            next_genotype += 1

    print(f"Cluster to genotype mapping: {cluster_to_genotype}")

    # 更新基因型
    for genotype in genotypes:
        if genotype == '.' or genotype == './.' or genotype == '.|.':
            new_genotypes.append(genotype)
        else:
            # 解析基因型
            is_phased = '|' in genotype
            separator = '|' if is_phased else '/'
            
            alleles = genotype.split(separator)
            new_alleles = []
            
            for allele in alleles:
                if allele == '.':
                    new_alleles.append('.')
                else:
                    try:
                        original_allele = int(allele)
                        if original_allele in cluster_mapping:
                            cluster_id = cluster_mapping[original_allele]
                            new_allele = str(cluster_to_genotype[cluster_id])
                            new_alleles.append(new_allele)
                            
                            # 统计AC（不包括REF）
                            if cluster_id != ref_cluster_id:
                                alt_idx = int(new_allele)
                                ac_counts[alt_idx] += 1
                        else:
                            # ALT没有映射，保持原样
                            new_alleles.append(allele)
                    except ValueError:
                        new_alleles.append(allele)
            
            new_genotype = separator.join(new_alleles)
            new_genotypes.append(new_genotype)

    print(f"Original genotypes: {genotypes}")
    print(f"New genotypes: {new_genotypes}")
    print(f"AC counts: {ac_counts}")

    # 更新ALT序列：选择每个簇的代表序列
    new_alts = []
    cluster_to_alt_seqs = {}
    
    # 收集每个簇的ALT序列（排除REF簇）
    for alt_id, cluster_id in cluster_mapping.items():
        if cluster_id != ref_cluster_id and alt_id > 0:
            if alt_id - 1 < len(alts):
                alt_seq = alts[alt_id - 1]
                if cluster_id not in cluster_to_alt_seqs:
                    cluster_to_alt_seqs[cluster_id] = []
                cluster_to_alt_seqs[cluster_id].append(alt_seq)

    # 为每个簇选择代表序列（按基因型顺序）
    for genotype_idx in sorted(cluster_to_genotype.values()):
        if genotype_idx == 0:
            continue  # 跳过REF
        
        # 找到对应的簇ID
        cluster_id = None
        for cid, gid in cluster_to_genotype.items():
            if gid == genotype_idx and cid in cluster_to_alt_seqs:
                cluster_id = cid
                break
        
        if cluster_id and cluster_id in cluster_to_alt_seqs:
            # 选择最长的序列作为代表
            longest_alt = max(cluster_to_alt_seqs[cluster_id], key=len)
            new_alts.append(longest_alt)

    print(f"New ALTs: {[f'len={len(alt)}bp' for alt in new_alts]}")

    # 更新INFO字段
    ac_values = [str(ac_counts[i]) for i in range(1, len(new_alts) + 1)]
    info_fields = info.split(';')
    updated_info_fields = []
    
    for field in info_fields:
        if field.startswith('AC='):
            updated_info_fields.append('AC=' + ','.join(ac_values))
        else:
            updated_info_fields.append(field)

    # 计算AN
    an = 0
    for genotype in genotypes:
        if genotype != '.' and genotype != './.' and genotype != '.|.':
            if '/' in genotype:
                alleles = genotype.split('/')
            elif '|' in genotype:
                alleles = genotype.split('|')
            else:
                alleles = [genotype]
            
            for allele in alleles:
                if allele != '.':
                    an += 1

    # 更新AN
    an_updated = False
    for i, field in enumerate(updated_info_fields):
        if field.startswith('AN='):
            updated_info_fields[i] = f'AN={an}'
            an_updated = True
            break
    if not an_updated:
        updated_info_fields.append(f'AN={an}')

    # 计算AF
    af_values = []
    for i in range(1, len(new_alts) + 1):
        if an > 0 and i in ac_counts:
            af = ac_counts[i] / an
            af_values.append(f"{af:.6f}")
        else:
            af_values.append("0.000000")
    
    # 更新AF
    af_updated = False
    for i, field in enumerate(updated_info_fields):
        if field.startswith('AF='):
            updated_info_fields[i] = 'AF=' + ','.join(af_values)
            af_updated = True
            break
    if not af_updated:
        updated_info_fields.append('AF=' + ','.join(af_values))
    
    updated_info = ';'.join(updated_info_fields)
    
    print(f"Updated INFO: {updated_info}")
    print("=" * 50)
    
    return (chrom, pos, id, ref, new_alts, updated_info, format_field, new_genotypes)

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
    with open(output_file, 'a') as file:
        for (chrom, pos), alt_dict in sorted(cluster_mapping.items()):
            for alt_id, cluster_id in alt_dict.items():
                file.write(f"{chrom}\t{pos}\t{alt_id}\t{cluster_id}\n")

def main(input_vcf, output_vcf, tmp_fasta, new_cluster_file, threads, 
         length_threshold=50, eps_short=0.5, eps_long=0.3, batch_size=1000):
    """主函数"""
    
    try:
        import numpy as np
        from sklearn.cluster import DBSCAN
    except ImportError as e:
        print(f"Error: Required package not found: {e}")
        print("Please install: pip install numpy scikit-learn")
        exit(1)
    
    # 读取VCF
    variants, header = read_vcf(input_vcf)
    write_vcf_header(output_vcf, header)
    
    print(f"Processing {len(variants)} variants")
    print(f"Clustering parameters: length_threshold={length_threshold}bp, eps_short={eps_short}, eps_long={eps_long}")
    
    # 处理每个变异
    for idx, variant in enumerate(variants):
        chrom, pos, id, ref, alts, info, format_field, genotypes = variant
        
        # 进度显示
        if idx % batch_size == 0:
            print(f"\nProcessed {idx}/{len(variants)} variants...")
        
        # 如果没有ALT或只有一个ALT，直接输出
        if len(alts) <= 1:
            write_vcf_variant(output_vcf, variant)
            continue
        
        # 写入临时FASTA
        write_tmp_fasta(chrom, pos, ref, alts, tmp_fasta)
        
        try:
            # 聚类（按长度分组）
            cluster_mapping, _ = cluster_sequences_with_length_groups(
                tmp_fasta, 
                length_threshold=length_threshold,
                eps_short=eps_short,
                eps_long=eps_long
            )
            
            # 保存聚类映射
            save_cluster_mapping({(chrom, pos): cluster_mapping}, new_cluster_file)
            
            # 更新变异信息，确保REF基因型为0
            updated_variant = update_genotypes_and_ac_with_ref_zero(variant, {(chrom, pos): cluster_mapping})
            write_vcf_variant(output_vcf, updated_variant)
            
        except Exception as e:
            print(f"\nError clustering {chrom}:{pos}: {e}")
            import traceback
            traceback.print_exc()
            write_vcf_variant(output_vcf, variant)
        
        finally:
            # 清理临时文件
            if os.path.exists(tmp_fasta):
                os.remove(tmp_fasta)
    
    print(f"\nProcessing completed. Output written to {output_vcf}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Cluster SV variants by length groups")
    parser.add_argument("input_vcf", help="Input VCF file")
    parser.add_argument("output_vcf", help="Output VCF file")
    parser.add_argument("tmp_fasta", help="Temporary FASTA file")
    parser.add_argument("new_cluster_file", help="File to store cluster mappings")
    parser.add_argument("--threads", type=int, default=1, help="Threads")
    parser.add_argument("--length_threshold", type=int, default=50, 
                       help="Length threshold for grouping (default: 50bp)")
    parser.add_argument("--eps_short", type=float, default=0.5, 
                       help="DBSCAN epsilon for short sequences (<threshold) (default: 0.5)")
    parser.add_argument("--eps_long", type=float, default=0.3, 
                       help="DBSCAN epsilon for long sequences (>=threshold) (default: 0.3)")
    parser.add_argument("--batch_size", type=int, default=1000, help="Progress interval (default: 1000)")
    
    args = parser.parse_args()
    
    main(args.input_vcf, args.output_vcf, args.tmp_fasta, args.new_cluster_file, 
         args.threads, args.length_threshold, args.eps_short, args.eps_long, args.batch_size)
