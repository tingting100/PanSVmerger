# 01_subsample_svs.sh
vcftools --gzvcf MC.pangraph.vcf.gz --indv sampleN --recode --stdout |bgzip -c > sampleN.sv.vcf.gz
tabix -p vcf sampleN.sv.vcf.gz

# 02_filter_large_svs.sh
vcfbub -l 0 -a 100000 -i sampleN.vcf.gz -o sampleN.filtered.vcf.gz

# 03_decompose_alleles.sh
vcfwave -I 1000 -i sampleN.filtered.vcf.gz -o sampleN.wave.vcf.gz

# 04_split_biallelic.sh
bcftools norm -m -any -i sampleN.wave.vcf.gz -o sampleN.biallelic.vcf.gz

# 05_run_truvari.sh
truvari bench -b sampleN.truth.sv.vcf.gz -c sampleN.biallelic.vcf.gz -o truvari_results
