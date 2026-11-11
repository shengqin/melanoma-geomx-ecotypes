"""Build Supplementary_Tables.xlsx (Supplementary Tables 1-4) from existing result tables.

Nothing is recomputed here; each sheet is a tidy copy of a table that already exists:
  Table S1  scRNA-seq marker genes            <- supplementary/source/SuppTable_scRNA_markers_all.csv
  Table S2  cell-state statistics by subtype  <- analysis/2026-11-10_state_by_subtype_all_AOIs/
  Table S3  per-AOI ecotypes and cell states  <- analysis/2026-01-28_GeoMx_recovery_noNorm/
  Table S4  cell-type-intrinsic genes         <- analysis/2026-06-02_hires_subtype_genes/scan_*.csv

Run:
uv run --with pandas --with openpyxl python scripts/build_supplementary_tables.py
"""
import pandas as pd

# Step 1. Paths.
BASE = "/Users/Brian/Library/CloudStorage/Box-Box/NanoString GeoMx data analysis"
FIGURES = f"{BASE}/manuscript/2026-09-15_manuscript/Figures"
RECOVERY = f"{BASE}/analysis/2026-01-28_GeoMx_recovery_noNorm"
STATE_STATS = f"{BASE}/analysis/2026-11-10_state_by_subtype_all_AOIs/state_by_subtype_all_vs_assigned.csv"
GENE_SCANS = f"{BASE}/analysis/2026-06-02_hires_subtype_genes"
OUTPUT = f"{FIGURES}/supplementary/Supplementary_Tables.xlsx"

# Display names used in the manuscript.
CELL_TYPE_NAME = {"Melanoma": "Melanoma", "T.CD8": "T-CD8", "T-CD8": "T-CD8", "T.CD4": "T-CD4", "T-CD4": "T-CD4",
                  "NK": "NK", "B_cell": "B cell", "Macrophage": "Macrophage", "Endothelial": "Endothelial",
                  "CAF": "CAF", "Mural_cell": "Mural cell"}
COMPARTMENT_NAME = {"S100": "S100+", "CD45": "CD45+", "neg/neg": "S100-/CD45-", "neg_neg": "S100-/CD45-"}
SUBTYPES = ["LMM", "SITU", "SSM", "NOD", "AMEL", "MET"]

# Step 2. Table S1: marker genes with adjusted P < 0.05 (the cut-off used for annotation),
# ranked within each population by adjusted P, then by log2 fold change.
markers = pd.read_csv(f"{FIGURES}/supplementary/source/SuppTable_scRNA_markers_all.csv")
markers = markers[markers["p_val_adj"] < 0.05].copy()
markers["population"] = markers["cluster"].map(CELL_TYPE_NAME)
markers = markers.sort_values(["population", "p_val_adj", "avg_log2FC"], ascending=[True, True, False])
table_s1 = markers[["population", "gene", "avg_log2FC", "pct.1", "pct.2", "p_val", "p_val_adj"]]
table_s1 = table_s1.rename(columns={"pct.1": "fraction_expressing_in_population",
                                    "pct.2": "fraction_expressing_in_other_cells",
                                    "p_val": "P_Wilcoxon", "p_val_adj": "P_adjusted_Bonferroni"})

# Step 3. Table S2: cell-state statistics. The all-AOI rows are the analysis reported in the paper;
# the assigned-only rows are kept for comparison.
state_stats = pd.read_csv(STATE_STATS)
state_stats["cell_type"] = state_stats["cell_type"].map(CELL_TYPE_NAME)
state_stats["compartment"] = state_stats["compartment"].map(COMPARTMENT_NAME)
state_stats["aoi_set"] = state_stats["aoi_set"].replace({"all_AOIs": "all melanoma AOIs (reported)",
                                                         "assigned_only": "ecotype-assigned AOIs only"})
ordered_columns = ["aoi_set", "compartment", "cell_type", "state", "n_aois", "n_patients"]
for subtype in SUBTYPES:
    ordered_columns.append(f"n_{subtype}")
    ordered_columns.append(f"patients_{subtype}")
    ordered_columns.append(f"median_{subtype}")
ordered_columns = ordered_columns + ["kruskal_p", "kruskal_q", "median_radial", "median_vertical",
                                     "radial_vs_vertical_p", "radial_vs_vertical_q",
                                     "patient_level_kruskal_p", "patient_level_kruskal_q",
                                     "patient_level_radial_vs_vertical_p", "patient_level_radial_vs_vertical_q"]
table_s2 = state_stats[ordered_columns]

# Step 4. Table S3: one row per AOI in the EcoTyper recovery input (296 nevus and melanoma AOIs).
metadata = pd.read_csv(f"{BASE}/data/2025/2025-01-27_GeoMxmetadata.txt", sep="\t")
metadata["ID"] = metadata["ID"].astype(str)
metadata = metadata.set_index("ID")
ecotype_abundance = pd.read_csv(f"{RECOVERY}/Ecotypes/ecotype_abundance.txt", sep="\t", index_col=0).T
ecotype_abundance.index = ecotype_abundance.index.astype(str)
assignment = pd.read_csv(f"{RECOVERY}/Ecotypes/ecotype_assignment.txt", sep="\t", index_col=0)
assignment.index = assignment.index.astype(str)
recovery_aois = list(ecotype_abundance.index)
table_s3 = pd.DataFrame(index=recovery_aois)
table_s3["slide"] = metadata.loc[recovery_aois, "Scan_ID"]
table_s3["ROI"] = metadata.loc[recovery_aois, "ROI..Label."]
table_s3["segment"] = metadata.loc[recovery_aois, "Segment_Tags"].map(COMPARTMENT_NAME)
table_s3["pathophysiology"] = metadata.loc[recovery_aois, "Pathophysiology"]
table_s3["subtype"] = metadata.loc[recovery_aois, "Path_sub_type"]
table_s3["patient"] = metadata.loc[recovery_aois, "Path_sub_type_patient"]
table_s3["ecotype_assignment"] = assignment["Ecotype"].reindex(recovery_aois).fillna("unassigned")
for ecotype in ecotype_abundance.columns:
    table_s3[f"{ecotype}_abundance"] = ecotype_abundance[ecotype]
for folder_name in ["Melanoma", "T.CD8", "T.CD4", "NK", "B_cell", "Macrophage", "Endothelial", "CAF", "Mural_cell"]:
    states = pd.read_csv(f"{RECOVERY}/{folder_name}/state_abundances.txt", sep="\t", index_col=0)
    states.columns = states.columns.astype(str)
    for state in states.index:
        column_name = f"{CELL_TYPE_NAME[folder_name]}-{state}"
        table_s3[column_name] = states.loc[state].reindex(recovery_aois)
table_s3.index.name = "AOI"
table_s3 = table_s3.reset_index()

# Step 5. Table S4: every gene tested per cell type, most significant first.
compartment_of = {"Melanoma": "S100", "CAF": "neg_neg", "Mural_cell": "neg_neg", "Endothelial": "neg_neg",
                  "T-CD4": "CD45", "T-CD8": "CD45", "NK": "CD45", "B_cell": "CD45", "Macrophage": "CD45"}
scans = []
for cell_type, compartment in compartment_of.items():
    scan = pd.read_csv(f"{GENE_SCANS}/scan_{cell_type}_{compartment}.csv")
    scan["cell_type"] = CELL_TYPE_NAME[cell_type]
    scan["compartment"] = COMPARTMENT_NAME[compartment]
    scan = scan.sort_values("p")
    scans.append(scan)
all_scans = pd.concat(scans, ignore_index=True)
table_s4 = all_scans[["cell_type", "compartment", "gene", "p", "q", "peak_subtype", "low_subtype", "range_log2", "epidermal"]]
table_s4 = table_s4.rename(columns={"p": "kruskal_p", "q": "q_BH_within_cell_type", "low_subtype": "lowest_subtype",
                                    "range_log2": "range_of_subtype_medians_log2", "epidermal": "epidermal_gene_excluded"})

# Step 6. A first sheet that says what each table is.
readme = pd.DataFrame({
    "Table": ["Supplementary Table 1", "Supplementary Table 2", "Supplementary Table 3", "Supplementary Table 4"],
    "Content": [
        "Marker genes of the nine populations of the integrated melanoma scRNA-seq reference (Seurat FindAllMarkers, "
        "Wilcoxon rank-sum test, positive markers, min.pct 0.25, log2FC threshold 0.25; adjusted P < 0.05).",
        "Cell-state abundance across melanoma subtypes in each GeoMx compartment. Kruskal-Wallis across six subtypes and "
        "Mann-Whitney radial (LMM, SITU, SSM) versus vertical (NOD, AMEL, MET); Benjamini-Hochberg q within compartment; "
        "patient-level = one median value per patient. Rows for ecotype-assigned AOIs only are given for comparison.",
        "Per-AOI annotation, EcoTyper ecotype assignment, ecotype abundance and cell-state abundances for the 296 nevus "
        "and melanoma AOIs.",
        "Cell-type-intrinsic genes (CIBERSORTx high-resolution expression) tested across the six melanoma subtypes within "
        "the compartment of each cell type; Kruskal-Wallis P and Benjamini-Hochberg q.",
    ],
    "Rows": [len(table_s1), len(table_s2), len(table_s3), len(table_s4)],
})

# Step 7. Write the workbook.
with pd.ExcelWriter(OUTPUT, engine="openpyxl") as writer:
    readme.to_excel(writer, sheet_name="README", index=False)
    table_s1.to_excel(writer, sheet_name="Table S1 scRNA markers", index=False)
    table_s2.to_excel(writer, sheet_name="Table S2 state statistics", index=False)
    table_s3.to_excel(writer, sheet_name="Table S3 per-AOI states", index=False)
    table_s4.to_excel(writer, sheet_name="Table S4 intrinsic genes", index=False)
print(readme[["Table", "Rows"]].to_string(index=False))
print("assigned AOIs in Table S3:", int((table_s3["ecotype_assignment"] != "unassigned").sum()), "of", len(table_s3))
print("genes with q < 0.05 and not epidermal, per cell type:")
significant = table_s4[(table_s4["q_BH_within_cell_type"] < 0.05) & (~table_s4["epidermal_gene_excluded"].astype(bool))]
print(significant["cell_type"].value_counts().to_string())
