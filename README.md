# SV-Merger-A-Pipeline-for-Merging-Multi-Allele-Structural-Variants
## About

SV-Merger is a bioinformatics pipeline designed to merge multi-allele structural variants (SVs) located at the same genomic position based on their length and sequence similarity. This tool aims to provide a reliable and efficient method for consolidating SV data, which is crucial for accurate genomic analysis and interpretation

## Requirements and Install



## Parameters



## Usage

examples

1. Xxx
2. xx
3. Xxx

## Output



## Citation



Xxx



## Help



If you get any problems, please raise an [issue](https://github.com/yeeus/GCI/issues) first. For another helps, please contact xxx



## TO-DO-LIST

- consider make most of parameters accessible, such as 
  - length threshold for grouping or make is as an dynamic setting which depends on allele length; 
  - the vsearch path, the main parameters of vsearch, like threads, cluster mode, 
  - for clustering very long or large sequences, consider to use [HAlign](https://github.com/malabz/HAlign-3)
  - 

- faster IO? using other method to read or write VCF files? [PyVCF](https://github.com/jamescasbon/PyVCF)
  - I will use `pysam` to edit vcf records more efficiently. (quanyu)

- Add some modules for statistics or visualization
  - basic statistics
  - histogram for AC before and after clustering

- other problem
  - whether should consider the direction of alignment, for example, in a clustered group, one sequence aligns with inverted direction, how to deal with it?
  - 
