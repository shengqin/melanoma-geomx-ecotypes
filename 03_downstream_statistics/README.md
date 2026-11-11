# 03. Downstream statistics

Python scripts that test the EcoTyper and CIBERSORTx outputs across the six melanoma subtypes
(LMM, SITU, SSM = radial growth; NOD, AMEL, MET = vertical growth or metastasis).

## `2026-11-10_state_by_subtype_all_AOIs.py` (Fig. 7d–f, Supplementary Table 2)
- **Input:** GeoMx metadata, EcoTyper recovery output (`<cell type>/state_abundances.txt`,
  `Ecotypes/ecotype_assignment.txt`) and the Q3-normalized GeoMx expression matrix.
- **Tests:** for every cell state in its own compartment (melanoma states in S100+ AOIs; T, NK, B and macrophage
  states in CD45+ AOIs; CAF, mural and endothelial states in S100−/CD45− AOIs):
  - Kruskal–Wallis across subtypes (subtypes with fewer than 3 AOIs left out);
  - Mann–Whitney, radial versus vertical growth;
  - Benjamini–Hochberg FDR across the states of each compartment;
  - the same two tests with one value per patient (median of that patient's AOIs), as a sensitivity check.
- **AOI sets:** all melanoma AOIs (reported in the paper) and, for comparison, only the AOIs EcoTyper assigned to an ecotype.
- **Also:** Spearman correlation of Melanoma-S01 with the pigmentation/MITF and IFN-γ signature scores in S100+ AOIs.
- **Output:** `state_by_subtype_all_vs_assigned.csv`, `melanoma_S01_vs_signatures.csv`.
- **Run:** `uv run --with pandas --with scipy --with numpy --with statsmodels python 2026-11-10_state_by_subtype_all_AOIs.py`

## `hires_all.py` (Fig. 8, Supplementary Table 4)
- **Input:** CIBERSORTx high-resolution expression for the GeoMx AOIs (`CIBERSORTxHiRes_NA_<cell type>_Window36.txt`) and the GeoMx metadata.
- **What it does:** for each cell type, keeps the melanoma AOIs of the compartment where that cell type is imputed
  (melanoma: S100+; CAF, endothelial, mural: S100−/CD45−; T-CD8, T-CD4, NK, B cell, macrophage: CD45+),
  drops genes not estimated in at least 40% of those AOIs, log2-transforms, and tests every gene across subtypes
  (Kruskal–Wallis, Benjamini–Hochberg FDR within the cell type). Epidermal and keratinocyte genes are flagged
  (`epidermal` column) so they can be excluded as spillover.
- **Output:** `scan_<cell type>_<compartment>.csv` (one per cell type), `candidate_genes_all.csv`, `summary_all.pkl`.

## Before running
Set the base path at the top of each script (`BASE`, and `OUT` in `hires_all.py`) to your own folders.
