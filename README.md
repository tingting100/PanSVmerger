# PanSVmerger: A Pipeline for Consolidating Multi-Allele Structural Variants

## 📖 Introduction
**PanSVmerger** is a high-performance bioinformatics pipeline designed to refine and merge structural variant (SV) calls derived from pangenome graphs (e.g., Minigraph-Cactus). 

It performs three main stages:
1.Regional filtering​ (optional, restrict to regions of interest)
2.SV locus extraction​ (select loci where any ALT differs from REF by ≥50 bp)
3.Intra-locus allele clustering via three strategies:
  Strategy A: k-mer-based adaptive Jaccard clustering
  Strategy B: VSEARCH global-alignment clustering
  Strategy C: length-based grouping

---

## 🛠 Prerequisites

### Software Dependencies
- [vcflib >= 1.0.3](https://github.com/vcflib/vcflib) (for `vcfwave` and `vcfcreatemulti`)
- [bcftools >= 1.18](https://github.com/samtools/bcftools)
- [vsearch >= 2.22.1](https://github.com/torognes/vsearch) (required for Strategy B)

### Python Environment
- Python 3.8+
- `pysam`, `scikit-learn`, `biopython`, `numpy`, `pandas`, `tqdm`

Install dependencies via pip:
```bash
pip install pysam scikit-learn biopython numpy pandas tqdm
```

---

## 🚀 Pipeline Workflow

### Step 1: Complex Region Masking
Biological "noise" in regions like centromeres, telomeres, and repetitive elements (STRs/VNTRs) can lead to unreliable SV calls. We filter these out using a BED-based approach.

```bash
python bin/filter_tr_region.py \
    --vcf standardized.vcf \
    --bed complex_regions.bed \
    --out filtered.vcf \
```

### Step 2: SV locus extraction
For VCF record, PanSVmerger identifies loci where the difference length between reference and alternate allele sequences is ≥50 bp.

```bash
perl bin/graph.getSV.pl \
    --vcf filtered.vcf \ 
    --out filtered.SV.vcf \
```

### Step 3: Multi-Allele Clustering
PanSVmerger offers three complementary clustering strategies to merge redundant alleles at the same locus into representative consensus alleles.

#### Strategy A: K-mer-based adaptive clustering (Jaccard method)

**Script:** `SV.weighted_joint_clustering.py`  
This strategy utilizes a hybrid distance matrix combining **K-mer frequency (Jaccard Distance)** and **Normalized Length Differences**. It uses DBSCAN to identify clusters and select a representative allele for each.
*   **Use Case:** Default choice for general pangenome analysis; balances sequence composition and variant size.

#### Strategy B: Global alignment-based clustering (Vsearch method)
**Script:** `SV.sequence.similarity.py`  
Leverages **VSEARCH** to perform global alignment between all alleles at a locus. Alleles are merged if they meet a specific identity threshold (e.g., 90%).
*   **Use Case:** When high-resolution nucleotide identity is critical.

#### Strategy C: Length distribution-based clustering (Length method)
**Script:** `SV.length.similarity.py`  
A rapid approach focusing solely on the **Length Distribution** of alleles using DBSCAN.
*   **Use Case:** Rapid processing of massive datasets where sequence composition is less important than SV size.

---

## 📊 Strategy Comparison

| Feature | Strategy A | Strategy B | Strategy C |
| :--- | :--- | :--- | :--- |
| **Primary Metric** | K-mer Count + $\Delta$ Length | Global Sequence % Identity | $\Delta$ Length (1D) |
| **Algorithm** | DBSCAN | VSEARCH | DBSCAN |
| **Computational Cost**| Medium | High | Low |
| **Clustering Logic** | Compositional Similarity | Base-level Precision | Size Distribution |

---

## ⚙️ Key Parameters

- `--eps`: (DBSCAN) The maximum distance between two alleles to be clustered together. A larger `eps` results in more aggressive merging.
- `--lambda`: (Weighted Strategy) A scaling factor that balances the weight between k-mer similarity and length difference.
- `--k`: (Weighted Strategy) K-mer size for Jaccard distance calculation (default: 5).
- `--padding`: (Filtering) Boundary extension (bp) around the masked BED regions.

---

## 📂 Project Structure
```text
PanSVMerger/
├── bin/
│   ├── filter_tr_region.py           # Masking complex regions
│   ├── SV.length.similarity.py        # Strategy C
│   ├── SV.sequence.similarity.py      # Strategy B
│   └── SV.weighted_joint_clustering.py # Strategy A
├── example/                          # Test datasets
├── README.md
└── requirements.txt
```

---

## 🤝 Contact & Citation
If you use PanSVmerger in your research, please cite:
PanSVmerger: a flexible pipeline for merging multiallelic structural variants in pangenome graphs

**Support:**
For questions or bug reports, please open an [issue](https://github.com/tingting100/PanSVMerger/issues) or contact:
- yangting@caas.cn
```
