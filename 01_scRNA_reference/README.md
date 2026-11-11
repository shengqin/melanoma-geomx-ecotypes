# 01. Single-cell melanoma reference

`2026-01-15_scRNA.R` builds the nine-population scRNA-seq reference and writes the inputs for CIBERSORTx
(Fig. 5a–c). It is an interactive analysis script: run it section by section in RStudio, not top to bottom.

## Input
- `2025-10-14_scRNA_melanoma.rds`: Seurat object of the three merged and Harmony-integrated melanoma datasets
  (GSE215121, GSE115978 and the 10x Genomics "Melanoma Tumor Derived Cells (v2)" dataset). The merge and
  integration that produced this object are not part of this script.

## What each section does
| Section | What it does |
|---|---|
| Start, "Subcluster for melanoma cell", "Annotate T-NK / B-Macro / fibro-endo" | Louvain clustering (mostly resolution 0.2) of the whole object, then sub-clustering of each lineage to annotate it. |
| "Filter", "Reannotate T-NK", "3rd time clean" | Removes low-quality and mixed clusters and re-clusters, in several rounds. |
| "Refine-melanoma cluster", "Refine-t-nk cluster" | Final clean-up of the melanoma and T/NK compartments. |
| Saved object | `2026-01-25_melanoma_filter.rds`: the final reference (48,733 cells, nine populations in `mega_cluster`). |
| "Make Pseudobulk and expression matrix" | Pseudobulk mixtures for the benchmark (Fig. 5c), one file per population: that population is added to a fixed background of all other cells at 0–40% in 2% steps (21 mixtures), averaging 10 samplings with replacement (seeds 1–10); counts are converted to CPM. Also writes `2026-01-28_melanoma_obj_feature_4000.txt`, the single-cell reference for the CIBERSORTx signature matrix (4,000 variable protein-coding genes). |
| "Analyze pseudo bulk result" | Reads the CIBERSORTx results for the mixtures and computes estimate-versus-truth correlations. |
| "Make fraction plot" | GeoMx melanoma fraction by segment (precursor of Fig. 5d). |
| "Make melanoma only data" | Removes normal-skin AOIs from the GeoMx matrix, giving the EcoTyper recovery input `2025-01-24_GeoMx_expression_data_noNorm.txt`. |
| Last sections | Exploratory heat map and PCA of melanoma high-resolution expression; not used in the paper. |

## Outputs used downstream
- `2026-01-28_melanoma_obj_feature_4000.txt` → `02_cibersortx_ecotyper_cluster_jobs/` (signature matrix and HiRes).
- `2026-01-28_pseudobulk_melanoma_<cell type>_cpm.txt` → the pseudobulk benchmark job.
- `2025-01-24_GeoMx_expression_data_noNorm.txt` → EcoTyper recovery.
- `2026-01-25_melanoma_filter.rds` → Fig. 5a and Supplementary Fig. 1 (`04_figures/fig5_supp_R_panels.R`).

## Before running
- Set `file_path` on line 6 to your data folder. A few later sections use relative paths (`../data/`).
- Packages: Seurat 5, harmony, dplyr, ggplot2, ggpubr, readr, DescTools, tidyr, limma, pheatmap.
