import os
import argparse
from collections import Counter
import shutil
import gzip
import re
## allele length ==1 not included: vsearch_path} --cluster_fast {input_fasta} --id {identity} --strand both --uc {output_prefix}.uc --threads 4 --maxseqlength 100000 --minseqlength 2

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
    #with gzip.open(vcf_file, 'rt') as file:
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

def run_vsearch(input_fasta, output_prefix, threads, identity=0.9):
    vsearch_path = "/share/home/zhanglab/user/yangting/software/vsearch/bin/vsearch"
    command = f"{vsearch_path} --cluster_fast {input_fasta} --id {identity}  --strand both --uc {output_prefix}.uc --threads {threads} --maxseqlength 100000 --minseqlength 2"
    print(f"Running command: {command}")
    os.system(command)


def parse_vsearch_clusters(uc_file):
    clusters = {}
    insufficient_fields_lines = []

    # Read the .uc file and parse clusters
    with open(uc_file, 'r') as file:
        for line in file:
            if line.startswith('S') or line.startswith('H'):
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    insufficient_fields_lines.append(line.strip())
                    continue

                try:
                    cluster_id = int(fields[1])
                    length = int(fields[2])
                    seq_id = fields[8]

                    # Check if this sequence is the reference
                    is_ref = '-ref' in seq_id

                    chrom, pos, alt_id = seq_id.split('-')
                    pos = int(pos)
                    alt_id = int(alt_id.replace("alt", "")) if "alt" in alt_id else 0

                    if (chrom, pos) not in clusters:
                        clusters[(chrom, pos)] = {}

                    if cluster_id not in clusters[(chrom, pos)]:
                        clusters[(chrom, pos)][cluster_id] = []

                    # Add the sequence info to the cluster, marking if it's a reference
                    clusters[(chrom, pos)][cluster_id].append((alt_id, length, is_ref))
                except (ValueError, IndexError) as e:
                    insufficient_fields_lines.append(line.strip())
                    continue

    # Step 1: Sort sequences in clusters by length and assign new sub-cluster IDs
    new_clusters = {}
    max_cluster_id = 0

    for (chrom, pos), cluster_dict in clusters.items():
        new_clusters[(chrom, pos)] = {}
        ref_cluster_id = None
        previous_length = None

        for cluster_id, seqs in cluster_dict.items():
            seqs.sort(key=lambda x: x[1], reverse=True)
            current_cluster_id = max_cluster_id
            for i, (alt_id, length, is_ref) in enumerate(seqs):
                if i == 0:
                    new_clusters[(chrom, pos)][alt_id] = current_cluster_id
                    previous_length = length
                    longest_alt_id = alt_id 
                else:
                    if previous_length - length > 50:
                        current_cluster_id += 1
                    new_clusters[(chrom, pos)][alt_id] = current_cluster_id
                    previous_length = length
                if is_ref:
                    ref_cluster_id = current_cluster_id
                    ref_alt_id = alt_id

            max_cluster_id = current_cluster_id + 1

        # Step 2: Set the reference cluster_id to 0 and adjust other cluster IDs
        if ref_cluster_id is not None:
            for alt_id in new_clusters[(chrom, pos)]:
                if new_clusters[(chrom, pos)][alt_id] == ref_cluster_id:
                    new_clusters[(chrom, pos)][alt_id] = 0
            
            cluster_id_mapping = {}
            new_cluster_id = 1
            for alt_id in sorted(new_clusters[(chrom, pos)]):
                old_cluster_id = new_clusters[(chrom, pos)][alt_id]
                if old_cluster_id != 0:
                    if old_cluster_id not in cluster_id_mapping:
                        cluster_id_mapping[old_cluster_id] = new_cluster_id
                        new_cluster_id += 1
                    new_clusters[(chrom, pos)][alt_id] = cluster_id_mapping[old_cluster_id]
        else:
            # 如果 ref_cluster_id 为空，则所有的 alt_id 的 cluster_id 都加1
            for alt_id in new_clusters[(chrom, pos)]:
                new_clusters[(chrom, pos)][alt_id] += 1
            new_clusters[(chrom, pos)][0] = 0  # 将 cluster_id 为0 的位置留给参考等位基因

    return new_clusters,insufficient_fields_lines

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

def update_genotypes_and_ac(variant, clusters):
    chrom, pos, id, ref, alts, info, format_field, genotypes = variant
    clusters = ensure_ref_and_alts_in_clusters(variant, clusters)
    cluster_mapping = clusters[(chrom, pos)]

    # Debug information: print cluster mapping
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

    # Update ALT sequences: take the longest sequence from each cluster
    new_alts = []
    cluster_to_alt = {cluster_id: [] for cluster_id in sorted(cluster_to_genotype.keys())}

    for alt_id, cluster_id in cluster_mapping.items():
        if cluster_id != ref_cluster_id:
            # Store alternative alleles based on their cluster ID
            cluster_to_alt[cluster_id].append(alts[alt_id - 1] if alt_id - 1 < len(alts) else ref)

    # Select the longest sequence for each alternative cluster
    for cluster_id in sorted(cluster_to_genotype.keys()):
        if cluster_id != ref_cluster_id:
            longest_alt = max(cluster_to_alt[cluster_id], key=len, default=ref)
            new_alts.append(longest_alt)
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
    # Update INFO field with new AN and NS value
    an_updated = False
    for i, field in enumerate(updated_info_fields):
        if field.startswith('AN='):
            updated_info_fields[i] = f'AN={an}'
            an_updated = True
        elif field.startswith('NS='):
            updated_info_fields[i] = f'NS={an}'
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

    # Update INFO field with new AT values
    for i, field in enumerate(updated_info_fields):
        if field.startswith('AT='):
            at_values = updated_info_fields[i].split(',')[:len(new_alts)]
            updated_info_fields[i] = ','.join(at_values)
            break
    
    updated_info = ';'.join(updated_info_fields)
    
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

def save_insufficient_fields_lines(insufficient_fields_lines, output_file):
    with open(output_file, 'w') as file:
        for line in insufficient_fields_lines:
            file.write(line + '\n')

def append_uc_to_merged(uc_file, merged_uc_file, lines_to_append, first_write=False):
    # 如果是第一次写入，清空文件内容
    if first_write:
        with open(merged_uc_file, 'w') as outfile:
            outfile.write('')

    # 以追加模式打开 merged_uc_file 并写入内容
    with open(merged_uc_file, 'a') as outfile:
        for line in lines_to_append:
            if line.strip():  # 避免追加空行
                outfile.write(line + '\n')

def save_cluster_mapping(cluster_mapping, output_file):
    # Save the cluster mapping to a file
    with open(output_file, 'a') as file:
        for (chrom, pos), alt_dict in sorted(cluster_mapping.items()):
            for alt_id, cluster_id in alt_dict.items():
                file.write(f"{chrom}\t{pos}\t{alt_id}\t{cluster_id}\n")

def main(input_vcf, output_vcf, tmp_fasta, cluster_prefix, merged_uc_file, new_cluster_file,insufficient_fields_file, threads):
    variants, header = read_vcf(input_vcf)
    write_vcf_header(output_vcf, header)
    first_write = True  # 标记是否是第一次写入 merged_uc_file
 
    for variant in variants:
        chrom, pos, id, ref, alts, info, format_field, genotypes = variant
        if len(alts) == 1:
            # 直接输出
            write_vcf_variant(output_vcf, variant)
        else:
            long_alts = [alt for alt in alts if len(alt) > 50000]
            if len(long_alts) > 100:
                # 直接输出原始变异信息
                write_vcf_variant(output_vcf, variant)
                continue

            write_tmp_fasta(chrom, pos, ref, alts, tmp_fasta)
            run_vsearch(tmp_fasta, cluster_prefix, threads)
            
            # Parse VSEARCH output to get clusters
            initial_clusters,insufficient_fields_lines = parse_vsearch_clusters(f"{cluster_prefix}.uc")
            # 如果没有聚类信息，直接输出原始变异信息
            if (chrom, pos) not in initial_clusters:
                write_vcf_variant(output_vcf, variant)
                continue

            # Ensure reference and alternative alleles are included in clusters
            updated_clusters = ensure_ref_and_alts_in_clusters(variant, initial_clusters)
            
            with open(f"{cluster_prefix}.uc", 'r') as infile:
                lines_to_append = infile.readlines()
            append_uc_to_merged(f"{cluster_prefix}.uc", merged_uc_file, lines_to_append, first_write=True)
            first_write = False
            save_cluster_mapping(updated_clusters, new_cluster_file)
            save_insufficient_fields_lines(insufficient_fields_lines, insufficient_fields_file)
            updated_variant = update_genotypes_and_ac(variant, updated_clusters)
            write_vcf_variant(output_vcf, updated_variant)
            
            # Clean up temporary files
#            for f in os.listdir():
 #               if f.startswith(cluster_prefix):
  #                  os.remove(f)
#

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process VCF file and cluster variants using VSEARCH.")
    parser.add_argument("input_vcf", help="Input VCF file")
    parser.add_argument("output_vcf", help="Output VCF file")
    parser.add_argument("tmp_fasta", help="Temporary FASTA file for VSEARCH")
    parser.add_argument("cluster_prefix", help="Prefix for VSEARCH cluster output files")
    parser.add_argument("merged_uc_file", help="File to store the merged cluster statistics")
    parser.add_argument("new_cluster_file", help="File to store the new cluster mappings")
    parser.add_argument("insufficient_fields_file", help="File to store not cluster mappings")
    parser.add_argument("threads", type=int, help="Threads to run vsearch")
    
    args = parser.parse_args()
    main(args.input_vcf, args.output_vcf, args.tmp_fasta, args.cluster_prefix, args.merged_uc_file, args.new_cluster_file,args.insufficient_fields_file, args.threads)
