# Melanoma GeoMx ecotypes: analysis code

Code for the single-cell reference, deconvolution, EcoTyper and figure analyses (Figs. 5–8 and the
supplementary figure and tables) of the melanoma GeoMx tissue-microarray study by Chhabra, Su et al.

No data are stored here. See "Data" below for where each input comes from.

## Folders

| Folder | What it does |
|---|---|
| `01_scRNA_reference/` | Integrates the three melanoma scRNA-seq datasets (Seurat + Harmony), annotates the nine populations, and writes the inputs for the CIBERSORTx signature matrix and the pseudobulk benchmark. |
| `02_cibersortx_ecotyper_cluster_jobs/` | SLURM jobs: CIBERSORTx signature matrix (400–800 genes per population), pseudobulk benchmark, high-resolution imputation, EcoTyper discovery on TCGA-SKCM and EcoTyper recovery on the GeoMx AOIs. |
| `03_downstream_statistics/` | Cell-state abundance across melanoma subtypes (all AOIs, with a patient-level check) and the cell-type-intrinsic gene scan. |
| `04_figures/` | Panels and page layout for Figs. 5–8 and Supplementary Fig. 1, and the supplementary tables workbook. |

## Before running

- The cluster jobs call CIBERSORTx, which needs an account. Replace `YOUR_CIBERSORTX_EMAIL` and
  `YOUR_CIBERSORTX_TOKEN` with your own credentials from https://cibersortx.stanford.edu.
- The scripts use absolute paths from the original analysis. Change the base path at the top of each script.

## Data

- scRNA-seq: GEO GSE215121 and GSE115978, and the 10x Genomics "Melanoma Tumor Derived Cells (v2)" dataset
  (https://www.10xgenomics.com/datasets/melanoma-tumor-derived-cells-v-2-2-standard-4-0-0).
- Bulk RNA-seq: TCGA-SKCM, from the Genomic Data Commons.
- GeoMx whole-transcriptome data: see the Data availability statement of the paper.

## Software

R with Seurat 5 and Harmony; CIBERSORTx (Fractions and HiRes); EcoTyper; Python 3 with pandas, numpy, scipy,
statsmodels, matplotlib and PyMuPDF.
