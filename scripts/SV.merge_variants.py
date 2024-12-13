import sys
import argparse
from collections import Counter
import pysam
import random
import string
import shutil
import subprocess
import os


def write_tmp_fasta(chrom, pos, ref, alts, tmp_fa):
    with open(tmp_fa, 'w') as file:
        file.write(f">{chrom}-{pos}-ref\n{ref}\n")
        for j, alt in enumerate(alts):
            file.write(f">{chrom}-{pos}-alt{j+1}\n{alt}\n")


def run_vsearch(input_fasta, threads, identity, maxseqlength, minseqlength):
    vsearch_path = shutil.which('vsearch')
    if not vsearch_path:
        sys.exit('ERROR! Please check `vesearch` is in your environment path\n')
    else:
        subprocess.run(f"{vsearch_path} --cluster_fast {input_fasta} --id {identity}  --strand both --uc {input_fasta}.uc --threads {threads} --maxseqlength {maxseqlength} --minseqlength {minseqlength}", shell=True, check=True)


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


def ensure_ref_and_alts_in_clusters(rec, clusters):
    # Initialize cluster dictionary for the variant's position if it doesn't exist
    if (rec.chrom, rec.pos) not in clusters:
        clusters[(rec.chrom, rec.pos)] = {}
    
    # Ensure the reference sequence (ID = 0) is in the cluster
    if 0 not in clusters[(rec.chrom, rec.pos)]:
        clusters[(rec.chrom, rec.pos)][0] = 0  # Add reference sequence with cluster ID 0

    # Add alternative alleles and ensure they have unique IDs
    for i in range(1, len(rec.alts) + 1):
        if i not in clusters[(rec.chrom, rec.pos)]:
            # Assign a new cluster ID ensuring uniqueness
            clusters[(rec.chrom, rec.pos)][i] = max(clusters[(rec.chrom, rec.pos)].values(), default=0) + 1
    
    return clusters


def append_uc_to_merged(merged_uc_file, lines_to_append, first_write=False):
    # 如果是第一次写入，清空文件内容
    if first_write:
        with open(merged_uc_file, 'w') as outfile:
            outfile.write('')

    # 以追加模式打开 merged_uc_file 并写入内容
    with open(merged_uc_file, 'a') as outfile:
        for line in lines_to_append:
            if line.strip():  # 避免追加空行
                outfile.write(line + '\n')


def save_cluster_mapping(cluster_mapping, output_file, first_write=False):
    if first_write:
        with open(output_file, 'w') as outfile:
            outfile.write('')

    # Save the cluster mapping to a file
    with open(output_file, 'a') as file:
        for (chrom, pos), alt_dict in sorted(cluster_mapping.items()):
            for alt_id, cluster_id in alt_dict.items():
                file.write(f"{chrom}\t{pos}\t{alt_id}\t{cluster_id}\n")


def save_insufficient_fields_lines(insufficient_fields_lines, output_file):
    with open(output_file, 'w') as file:
        for line in insufficient_fields_lines:
            if line.strip():  # 避免追加空行
                file.write(line + '\n')


def update_record(rec, clusters, vcf_header, vcf_out):
    #clusters = ensure_ref_and_alts_in_clusters(rec, clusters)
    cluster_mapping = clusters[(rec.chrom, rec.pos)]

    # Debug information: print cluster mapping
    #print(f"Cluster mapping for {chrom}-{pos}: {cluster_mapping}")

    new_genotypes = []
    ac_counts = Counter()

    # Find the reference sequence's cluster ID
    ref_cluster_id = cluster_mapping[0]
    
    # Redefine genotypes based on clusters
    cluster_to_genotype = {ref_cluster_id: 0}
    next_genotype = 1

    for cluster_id in sorted(set(cluster_mapping.values())):
        if cluster_id != ref_cluster_id:
            cluster_to_genotype[cluster_id] = next_genotype # int
            next_genotype += 1

    an_number = 0
    for value_dict in rec.samples.values():
        genotype = value_dict['GT'][0]
        if genotype == None:
            new_genotype = genotype
        else:
            if genotype in cluster_mapping:
                new_genotype = cluster_to_genotype[cluster_mapping[genotype]]
                if cluster_mapping[genotype] != ref_cluster_id: # int
                    ac_counts[new_genotype] += 1
            else:
                # Handle singleton sequences
                new_genotype = next_genotype
                next_genotype += 1
                ac_counts[new_genotype] += 1
            an_number += 1
        new_genotypes.append(new_genotype)

    # Update ALT sequences: take the longest sequence from each cluster
    new_alts = []
    cluster_to_alt = {cluster_id: [] for cluster_id in sorted(cluster_to_genotype.keys())}

    for alt_id, cluster_id in cluster_mapping.items():
        if cluster_id != ref_cluster_id:
            # Store alternative alleles based on their cluster ID
            cluster_to_alt[cluster_id].append(rec.alts[alt_id - 1] if alt_id - 1 < len(rec.alts) else rec.ref)

    # Select the longest sequence for each alternative cluster
    for cluster_id in sorted(cluster_to_genotype.keys()):
        if cluster_id != ref_cluster_id:
            longest_alt = max(cluster_to_alt[cluster_id], key=len, default=rec.ref)
            new_alts.append(longest_alt)
    # 如果 new_alts 为空，直接返回原始 variant
    if not new_alts:
        vcf_out.write(rec)
        return

    new_info = dict(rec.info)
    # Update INFO field with new AC, AF values
    new_info['AC'] = tuple([ac_counts[i] for i in range(1, len(new_alts) + 1)])
    new_info['AF'] = tuple([ac_counts[i] / new_info['AN'] for i in range(1, len(new_alts) + 1)])
    new_info['AT'] = new_info['AT'][:len(new_alts) + 1] ## simply extract the first `len(new_alts) + 1` allele traversals
    
    new_rec = vcf_header.new_record(contig=rec.chrom, start=rec.start, stop=rec.stop, alleles=[rec.ref] + new_alts, id=rec.id, qual=rec.qual, filter=rec.filter, info=new_info)
    vcf_out.write(new_rec)


def main(input_vcf, output_prefix, alt_length=50000, alt_number=100, identity=0.9, maxseqlength=100000, minseqlength=2, threads=1):
    if input_vcf.endswith('.vcf.gz'):
        suffix = 'vcf.gz'
    elif input_vcf.endswith('.vcf'):
        suffix = 'vcf'
    else:
        print(f'Warning! Input vcf has strange filename extension instead of ".vcf" or ".vcf.gz". Ignore it anyway.\n')
        suffix = 'vcf'
    vcf = pysam.VariantFile(input_vcf, threads=threads)
    vcf_out = pysam.VariantFile(f'{output_prefix}.{suffix}', 'w', header=vcf.header, threads=threads)
    
    for rec in vcf.fetch():
        if len(rec.alts) == 1:
            vcf_out.write(rec)
        else:
            long_alts = [alt for alt in rec.alts if len(alt) >= alt_length]
            if len(long_alts) > alt_number:
                vcf_out.write(rec)
                continue
            

            random_string = ''.join(random.choices(string.ascii_letters + string.digits, k=10))
            tmp_fasta = f'{output_prefix}.{random_string}.fa'
            write_tmp_fasta(rec.chrom, rec.pos, rec.ref, rec.alts, tmp_fasta)
            run_vsearch(tmp_fasta, threads, identity, maxseqlength, minseqlength)
            # Parse VSEARCH output to get clusters
            initial_clusters, insufficient_fields_lines = parse_vsearch_clusters(f"{tmp_fasta}.uc")
            # 如果没有聚类信息，直接输出原始变异信息
            if (rec.chrom, rec.pos) not in initial_clusters:
                vcf_out.write(rec)
                continue

            # Ensure reference and alternative alleles are included in clusters
            updated_clusters = ensure_ref_and_alts_in_clusters(rec, initial_clusters)

            # with open(f"{tmp_fasta}.uc", 'r') as infile:
            #     lines_to_append = infile.readlines()
            # append_uc_to_merged(f"{tmp_fasta}.uc.merge", lines_to_append, first_write=True)

            # save_cluster_mapping(updated_clusters, f"{tmp_fasta}.uc.merge.new", first_write=True)
            # save_insufficient_fields_lines(insufficient_fields_lines, f'{output_prefix}.uncluster.vcf')
            update_record(rec, updated_clusters, vcf.header, vcf_out)
            
            os.remove(tmp_fasta)
            os.remove(f"{tmp_fasta}.uc")
    vcf.close()
    vcf_out.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process VCF file from MC and cluster alleles using VSEARCH", formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("input_vcf", help="Input xxx.wave.bi.createmulti.vcf(.gz) file\nIf input compressed vcf file, please index first")
    parser.add_argument("output_prefix", help="Prefix of output files\nCompressed file in, compressed files out, and vice versa")
    parser.add_argument('-al', '--alt-length', metavar='INT', type=int, help='Skip records with alternative alleles >= INT bp [50000]', default=50000)
    parser.add_argument('-an', '--alt-number', metavar='INT', type=int, help='Skip records with > INT longer alternative alleles `--alt-length` [100]', default=100)
    parser.add_argument('-id', '--identity', metavar='FLOAT', type=float, help='Cluster alleles above FLOAT identity [0.9]', default=0.9)
    parser.add_argument('-maxl', '--maxseqlength', metavar='INT', type=int, help='Maximum sequence length for `vsearch` [100000]', default=100000)
    parser.add_argument('-minl', '--minseqlength', metavar='INT', type=int, help='Minimum sequence length for `vsearch` [2]', default=2)
    parser.add_argument('-t', '--threads', metavar='INT', type=int, help='Number of threads [1]', default=1)
    
    args = parser.parse_args()
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()
    main(args.input_vcf, args.output_prefix, args.alt_length, args.alt_number, args.identity, args.maxseqlength, args.minseqlength, args.threads)
