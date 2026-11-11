"""Cell-state abundance across melanoma subtypes: all AOIs versus EcoTyper-assigned AOIs only.

Why this script exists (10 Nov 2026):
The June analysis compared cell-state abundances across subtypes using only the AOIs that EcoTyper
assigned to an ecotype (128 of 198 melanoma AOIs). EcoTyper estimates state abundances for every AOI,
so the ecotype threshold is not needed for these tests. Gagan asked why the AOI numbers were small,
and the authors agreed to use all melanoma AOIs. This script runs the same tests on both AOI sets so
the old and new numbers can be compared side by side.

Tests (unchanged from 2026-06-02_immune_stromal_subtype_analysis.py):
- Kruskal-Wallis across the six subtypes (subtypes with fewer than 3 AOIs are left out).
- Mann-Whitney, radial growth (LMM, SITU, SSM) versus vertical growth (NOD, AMEL, MET).
- Benjamini-Hochberg FDR within each compartment, over every state of that compartment's cell types.

Run:
uv run --with pandas --with scipy --with numpy --with statsmodels python 2026-11-10_state_by_subtype_all_AOIs.py
"""
import os

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

# Step 1. Paths.
BASE = "/Users/Brian/Library/CloudStorage/Box-Box/NanoString GeoMx data analysis"
RECOVERY = f"{BASE}/analysis/2026-01-28_GeoMx_recovery_noNorm"
OUT = f"{BASE}/analysis/2026-11-10_state_by_subtype_all_AOIs"
os.makedirs(OUT, exist_ok=True)

# Step 2. Fixed orders and groupings.
# Subtype order runs from radial to vertical growth, as in Fig. 7.
SUBTYPES = ["LMM", "SITU", "SSM", "NOD", "AMEL", "MET"]
GROWTH_PHASE = {"LMM": "radial", "SITU": "radial", "SSM": "radial",
                "NOD": "vertical", "AMEL": "vertical", "MET": "vertical"}
# Each compartment is tested only for the cell types that belong to it (same families as June).
COMPARTMENT_CELL_TYPES = {
    "S100": ["Melanoma"],
    "CD45": ["T.CD8", "T.CD4", "Macrophage", "NK", "B_cell"],
    "neg/neg": ["CAF", "Mural_cell", "Endothelial"],
}
# A subtype needs at least this many AOIs to enter the Kruskal-Wallis test (June rule).
MIN_AOIS_PER_SUBTYPE = 3

# Step 3. GeoMx metadata for all 316 AOIs; keep the 198 melanoma AOIs.
metadata = pd.read_csv(f"{BASE}/data/2025/2025-01-27_GeoMxmetadata.txt", sep="\t")
metadata["ID"] = metadata["ID"].astype(str)
metadata = metadata.set_index("ID")
melanoma_metadata = metadata[metadata["Pathophysiology"] == "Melanoma"]

# Step 4. The 167 AOIs that EcoTyper assigned to an ecotype.
assignment = pd.read_csv(f"{RECOVERY}/Ecotypes/ecotype_assignment.txt", sep="\t", index_col=0)
assigned_ids = set(assignment.index.astype(str))

# Step 5. State abundances (AOIs x states) for every cell type. EcoTyper wrote these for all 296
# nevus and melanoma AOIs in the recovery input.
all_cell_types = []
for cell_types in COMPARTMENT_CELL_TYPES.values():
    for cell_type in cell_types:
        all_cell_types.append(cell_type)
state_abundance = {}
for cell_type in all_cell_types:
    table = pd.read_csv(f"{RECOVERY}/{cell_type}/state_abundances.txt", sep="\t", index_col=0)
    table.columns = table.columns.astype(str)
    state_abundance[cell_type] = table.T


def test_compartment(compartment, aoi_ids, aoi_set_name):
    """Run the subtype tests for every state of one compartment on the given AOIs."""
    rows = []
    for cell_type in COMPARTMENT_CELL_TYPES[compartment]:
        states = state_abundance[cell_type]
        for state in states.columns:
            # Step a. Abundance, subtype and patient of each AOI.
            table = pd.DataFrame({
                "value": states[state].reindex(aoi_ids),
                "subtype": melanoma_metadata.loc[aoi_ids, "Path_sub_type"],
                "patient": melanoma_metadata.loc[aoi_ids, "Path_sub_type_patient"],
            })
            table = table.dropna(subset=["value"])

            # Step b. Kruskal-Wallis across subtypes with enough AOIs.
            groups = []
            row = {"aoi_set": aoi_set_name, "compartment": compartment, "cell_type": cell_type, "state": state,
                   "n_aois": len(table), "n_patients": table["patient"].nunique()}
            for subtype in SUBTYPES:
                values = table.loc[table["subtype"] == subtype, "value"]
                row[f"n_{subtype}"] = len(values)
                row[f"patients_{subtype}"] = table.loc[table["subtype"] == subtype, "patient"].nunique()
                if len(values) >= MIN_AOIS_PER_SUBTYPE:
                    row[f"median_{subtype}"] = values.median()
                    groups.append(values.values)
                else:
                    row[f"median_{subtype}"] = np.nan
            row["kruskal_p"] = stats.kruskal(*groups).pvalue

            # Step c. Radial versus vertical growth (two-sided Mann-Whitney).
            phase = table["subtype"].map(GROWTH_PHASE)
            radial_values = table.loc[phase == "radial", "value"]
            vertical_values = table.loc[phase == "vertical", "value"]
            row["median_radial"] = radial_values.median()
            row["median_vertical"] = vertical_values.median()
            row["radial_vs_vertical_p"] = stats.mannwhitneyu(radial_values, vertical_values,
                                                             alternative="two-sided").pvalue

            # Step d. Patient-level sensitivity check: one value per patient (median of that patient's
            # AOIs), then the same two tests. AOIs from one patient are not independent, so this shows
            # whether a result depends on counting each AOI separately.
            per_patient = table.groupby("patient").agg(value=("value", "median"), subtype=("subtype", "first"))
            patient_groups = []
            for subtype in SUBTYPES:
                patient_values = per_patient.loc[per_patient["subtype"] == subtype, "value"]
                if len(patient_values) >= MIN_AOIS_PER_SUBTYPE:
                    patient_groups.append(patient_values.values)
            row["patient_level_kruskal_p"] = stats.kruskal(*patient_groups).pvalue
            patient_phase = per_patient["subtype"].map(GROWTH_PHASE)
            row["patient_level_radial_vs_vertical_p"] = stats.mannwhitneyu(
                per_patient.loc[patient_phase == "radial", "value"],
                per_patient.loc[patient_phase == "vertical", "value"], alternative="two-sided").pvalue
            rows.append(row)

    # Step e. Benjamini-Hochberg FDR across all states of this compartment.
    result = pd.DataFrame(rows)
    result["kruskal_q"] = multipletests(result["kruskal_p"], method="fdr_bh")[1]
    result["radial_vs_vertical_q"] = multipletests(result["radial_vs_vertical_p"], method="fdr_bh")[1]
    result["patient_level_kruskal_q"] = multipletests(result["patient_level_kruskal_p"], method="fdr_bh")[1]
    result["patient_level_radial_vs_vertical_q"] = multipletests(
        result["patient_level_radial_vs_vertical_p"], method="fdr_bh")[1]
    return result


# Step 6. Run both AOI sets for the three compartments.
all_results = []
for compartment in COMPARTMENT_CELL_TYPES:
    compartment_ids = melanoma_metadata.index[melanoma_metadata["Segment_Tags"] == compartment]
    assigned_compartment_ids = [aoi for aoi in compartment_ids if aoi in assigned_ids]
    all_results.append(test_compartment(compartment, list(compartment_ids), "all_AOIs"))
    all_results.append(test_compartment(compartment, assigned_compartment_ids, "assigned_only"))
results = pd.concat(all_results, ignore_index=True)
results.to_csv(f"{OUT}/state_by_subtype_all_vs_assigned.csv", index=False)

# Step 7. Print a compact comparison.
pd.set_option("display.width", 250)
pd.set_option("display.max_columns", 40)
shown_columns = ["aoi_set", "compartment", "cell_type", "state", "n_aois", "n_patients", "kruskal_p", "kruskal_q",
                 "median_radial", "median_vertical", "radial_vs_vertical_p", "radial_vs_vertical_q",
                 "patient_level_kruskal_p", "patient_level_radial_vs_vertical_p", "patient_level_radial_vs_vertical_q"]
print(results[shown_columns].to_string(index=False, float_format=lambda value: f"{value:.3g}"))
median_columns = ["aoi_set", "compartment", "cell_type", "state"]
for subtype in SUBTYPES:
    median_columns.append(f"median_{subtype}")
    median_columns.append(f"n_{subtype}")
    median_columns.append(f"patients_{subtype}")
print(results[median_columns].to_string(index=False, float_format=lambda value: f"{value:.2f}"))

# Step 8. Melanoma-S01 versus the pigmentation and IFN-gamma signatures in S100+ melanoma AOIs
# (Spearman), for both AOI sets. Signature scores are computed exactly as in June: mean z-score of
# log2(Q3 + 1) across the signature genes, z-scored over all AOIs in the expression file.
expression = pd.read_csv(f"{BASE}/data/2025/2025-01-24_GeoMx_expression_data.txt", sep="\t", index_col=0)
expression.columns = expression.columns.astype(str)
log_expression = np.log2(expression + 1)
z_scores = log_expression.sub(log_expression.mean(axis=1), axis=0).div(log_expression.std(axis=1) + 1e-9, axis=0)
PIGMENT_GENES = ["PMEL", "MLANA", "TYR", "TYRP1", "DCT", "MITF", "OCA2", "GPR143", "TRPM1", "SLC45A2", "MFSD12", "TSPAN10"]
IFN_GENES = ["STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11", "GBP1", "IDO1", "HLA-DRA", "CD8A", "GZMB", "PRF1", "CD274", "TAP1"]
signature_rows = []
s100_ids = list(melanoma_metadata.index[melanoma_metadata["Segment_Tags"] == "S100"])
assigned_s100_ids = [aoi for aoi in s100_ids if aoi in assigned_ids]
for aoi_set_name, aoi_ids in [("all_AOIs", s100_ids), ("assigned_only", assigned_s100_ids)]:
    for signature_name, genes in [("pigmentation_MITF", PIGMENT_GENES), ("IFN_gamma", IFN_GENES)]:
        genes_present = [gene for gene in genes if gene in z_scores.index]
        score = z_scores.loc[genes_present].mean(axis=0)
        paired = pd.DataFrame({"S01": state_abundance["Melanoma"]["S01"].reindex(aoi_ids),
                               "score": score.reindex(aoi_ids)}).dropna()
        rho, p_value = stats.spearmanr(paired["S01"], paired["score"])
        signature_rows.append({"aoi_set": aoi_set_name, "signature": signature_name, "genes_used": len(genes_present),
                               "n_aois": len(paired), "spearman_rho": rho, "p": p_value})
signature_results = pd.DataFrame(signature_rows)
signature_results.to_csv(f"{OUT}/melanoma_S01_vs_signatures.csv", index=False)
print(signature_results.to_string(index=False))
