# PanSVMerger: A Pipeline for Merging Multi-Allelic Structural Variants

## 📖 About
**PanSVMerger** 是一款专门为泛基因组图形（Pangenome Graphs，如 Minigraph-Cactus 输出）设计的结构变异（SV）精简与聚类工具。

在从图形中提取变异时，复杂区域往往会产生大量冗余、嵌套或极度相似的多等位基因（Multi-allelic sites）。PanSVMerger 通过“分解-过滤-聚类”的策略，将同一基因组位置上长度相似或序列一致的等位基因合并，生成高质量、无冗余的共识变异集，极大地简化了下游的群体遗传学分析。

---

## 🚀 Pipeline Overview

1.  **Decomposition**: 利用 `vcfwave` 分解复杂的图形气泡。
2.  **Normalization**: 使用 `vcfcreatemulti` 重建多等位基因记录。
3.  **Filtering**: 使用 `filter_tr_region.py` 过滤掉串联重复序列（TR）或其他复杂干扰区域。
4.  **Clustering**: 根据研究需求，从三种聚类策略中选择一种执行合并。

---

## 🛠 Requirements & Installation

### Dependencies
- **Tools**: [vcflib](https://github.com/vcflib/vcflib), [bcftools](https://github.com/samtools/bcftools), [vsearch](https://github.com/torognes/vsearch)
- **Python Libraries**: `pysam`, `scikit-learn`, `biopython`, `numpy`, `pandas`, `tqdm`

### Install
```bash
git clone https://github.com/tingting100/PanSVMerger.git
cd PanSVMerger
pip install -r requirements.txt 
