threads=60

date

source  ~/softwares/cactus-bin-v2.9.8/venv-cactus-v2.9.8/bin/activate

time cactus-pangenome ./js  HG002.list   \
--outDir ./MC_CHM13v2_HG002 \
--outName CHM13v2_HG002 \
--reference CHM13v2m \
--gbz \
--giraffe \
--vcf \
--chrom-vg \
--maxCores 60 \
--maxMemory 1000G \
--workDir ./MC_CHM13v2_HG002
