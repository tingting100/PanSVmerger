# PanSVmerger: A Pipeline for Consolidating Multi-Allele Structural Variants

## 📖 Introduction
**PanSVmerger** is a high-performance bioinformatics pipeline designed to refine and merge structural variant (SV) calls derived from pangenome graphs (e.g., Minigraph-Cactus). 

Structural variants extracted from graph bubbles often contain redundant, nested, or slightly divergent alleles due to sequencing errors, alignment artifacts, or graph complexity. PanSVmerger addresses these challenges by implementing a three-stage workflow: **Decomposition**, **Biological Region Filtering**, and **Intelligent Allele Clustering**. It aims to provide a clean, non-redundant consensus VCF suitable for downstream population genetics and comparative genomics.

---

## 🛠 Prerequisites

### Software Dependencies
- [vcflib](https://github.com/vcflib/vcflib) (for `vcfwave` and `vcfcreatemulti`)
- [bcftools](https://github.com/samtools/bcftools)
- [vsearch](https://github.com/torognes/vsearch) (required for Strategy B)

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
> *Xxx, et al. (2025). PanSVmerger: A pipeline for structural variant consolidation in pangenome graphs.*

**Support:**
For questions or bug reports, please open an [issue](https://github.com/tingting100/PanSVMerger/issues) or contact:
- yangting@caas.cn
```
