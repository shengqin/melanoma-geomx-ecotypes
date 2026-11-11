"""Panels for Figs. 6, 7 and 8, drawn at their exact printed size.

Typography standard (author request, 30 Oct 2026), same as Fig. 5:
  Arial; 7 pt for titles (panel, axis and legend titles);
  6 pt for tick labels and all other small text (legend entries, in-plot labels, statistics).
Every PDF is saved at the size it is printed, so compose.py places it unscaled.

Run:  python panels_fig6_8.py <output folder> [panel names ...]
Data sources are the same tables used for the September 2026 figures; the printed
statistics should match those in the manuscript.
Change on 10 Nov 2026: Figs. 7d-f use all melanoma AOIs (68 S100+, 68 CD45+, 62 stromal), not only
the AOIs EcoTyper assigned to an ecotype. Full statistics, including FDR and a patient-level check, are in
analysis/2026-11-10_state_by_subtype_all_AOIs/.
"""
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib._mathtext as mathtext_engine
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch, Rectangle
from scipy import stats

# Step 1. Paths.
BASE = "/Users/Brian/Library/CloudStorage/Box-Box/NanoString GeoMx data analysis"
ANALYSIS = f"{BASE}/analysis"
DISCOVERY = f"{ANALYSIS}/2026-01-28_TCGA_SKCM_Discovery_997_4_4"      # TCGA-SKCM EcoTyper discovery run
RECOVERY = f"{ANALYSIS}/2026-01-28_GeoMx_recovery_noNorm"             # GeoMx EcoTyper recovery run
HIRES = f"{BASE}/data/2026-01-28_melanoma_GeoMx_hires"                # CIBERSORTx high-resolution output
SCAN = f"{ANALYSIS}/2026-06-02_hires_subtype_genes"                   # per-cell-type subtype scans
OUTPUT_DIR = sys.argv[1] if len(sys.argv) > 1 else "panels"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Step 2. Typography. Math text (P-value exponents, italic gene names) is set in Arial, and the
# superscript shrink factor is set to 1 so exponents print at the same size as the surrounding
# text (otherwise matplotlib shrinks them to 70%, below the 6 pt minimum).
TITLE_PT = 7
SMALL_PT = 6
POINTS_PER_INCH = 72
mathtext_engine.SHRINK_FACTOR = 1.0
plt.rcParams.update({
    "font.family": "Arial",
    "font.size": SMALL_PT,
    "axes.titlesize": TITLE_PT,
    "axes.labelsize": TITLE_PT,
    "xtick.labelsize": SMALL_PT,
    "ytick.labelsize": SMALL_PT,
    "legend.fontsize": SMALL_PT,
    "legend.title_fontsize": TITLE_PT,
    "mathtext.fontset": "custom",
    "mathtext.rm": "Arial",
    "mathtext.it": "Arial:italic",
    "mathtext.bf": "Arial:bold",
    "mathtext.default": "regular",
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "xtick.major.size": 2.5,
    "ytick.major.size": 2.5,
    "pdf.fonttype": 42,          # embed TrueType so text stays editable
    "savefig.dpi": 300,
})

# Step 3. Shared names, orders and colours (same as Fig. 5).
CELL_TYPE_ORDER = ["Melanoma", "T-CD8", "T-CD4", "NK", "B_cell", "Macrophage", "Endothelial", "CAF", "Mural_cell"]
DISPLAY_NAME = {"Melanoma": "Melanoma", "T-CD8": "T-CD8", "T-CD4": "T-CD4", "NK": "NK", "B_cell": "B cell",
                "Macrophage": "Macrophage", "Endothelial": "Endothelial", "CAF": "CAF", "Mural_cell": "Mural cell",
                "T.CD8": "T-CD8", "T.CD4": "T-CD4"}
CELL_COLOUR = {"Melanoma": "#7B4A23", "T-CD8": "#2C6FB0", "T-CD4": "#7FB2E0", "NK": "#9C27B0", "B_cell": "#3FA34D",
               "Macrophage": "#E08A2B", "Endothelial": "#D7263D", "CAF": "#C85A54", "Mural_cell": "#6D597A"}
SUBTYPES = ["LMM", "SITU", "SSM", "NOD", "AMEL", "MET"]
VERTICAL_GROWTH = {"LMM": False, "SITU": False, "SSM": False, "NOD": True, "AMEL": True, "MET": True}
RADIAL_COLOUR = "#3E7C7C"
VERTICAL_COLOUR = "#C85A54"
E1_COLOUR = "#D7263D"
E2_COLOUR = "#F2C14E"
PRIMARY_COLOUR = "#F2C14E"
METASTASIS_COLOUR = "#8C8C8C"
MISSING_COLOUR = "#E6E6E6"


# ------------------------------------------------------------------ helpers


def vector_heatmap(axis, values, cmap, vmin=None, vmax=None):
    # Draw a matrix as vector rectangles (pcolormesh) on imshow's coordinate grid: cell (row, column) is
    # centred at (column, row) and row 0 is at the top. Edges are drawn in the cell's own colour (0.1 pt)
    # so PDF viewers show no hairline gaps between neighbouring cells.
    values = np.atleast_2d(values)
    n_rows, n_columns = values.shape
    column_edges = np.arange(n_columns + 1) - 0.5
    row_edges = np.arange(n_rows + 1) - 0.5
    mesh = axis.pcolormesh(column_edges, row_edges, values, cmap=cmap, vmin=vmin, vmax=vmax,
                           edgecolors="face", linewidth=0.1, rasterized=False)
    axis.set_xlim(-0.5, n_columns - 0.5)
    axis.set_ylim(n_rows - 0.5, -0.5)
    return mesh


def keep_colorbar_vector(colorbar):
    # matplotlib rasterizes colour bars with many colour steps; switch that off so the PDF stays vector.
    colorbar.solids.set_rasterized(False)
    return colorbar


def new_figure(width_pt, height_pt):
    # A figure exactly as large as the printed panel.
    return plt.figure(figsize=(width_pt / POINTS_PER_INCH, height_pt / POINTS_PER_INCH), layout="constrained")


def save_panel(figure, name):
    # Save the vector PDF (used by compose.py) and a PNG preview; never crop, so the size stays exact.
    pdf_path = os.path.join(OUTPUT_DIR, f"{name}.pdf")
    figure.savefig(pdf_path)
    figure.savefig(pdf_path.replace(".pdf", ".png"), dpi=150)
    plt.close(figure)
    print("saved", name)


def format_p(p_value):
    # P values: two significant digits; below 0.001 as a power of ten with a one-decimal mantissa.
    if p_value < 0.001:
        exponent = int(np.floor(np.log10(p_value)))
        mantissa = p_value / 10 ** exponent
        return f"{mantissa:.1f} × 10$^{{\\mathregular{{−{abs(exponent)}}}}}$"
    return f"{p_value:.2g}"


def italic(text):
    # Italic text (gene names) through math text, which is set in Arial Italic.
    return "$\\it{" + text + "}$"


def title_with_stats(axis, title, stats_lines):
    # 7 pt title above one or more 6 pt statistics lines, all centred over the axes.
    # The statistics sit 2 pt above the axes; the title is lifted by their height plus 2 pt.
    line_height_pt = 9          # 6 pt text plus room for a raised exponent
    stats_text = "\n".join(stats_lines)
    axis.annotate(stats_text, xy=(0.5, 1.0), xycoords="axes fraction", xytext=(0, 2), textcoords="offset points",
                  ha="center", va="bottom", fontsize=SMALL_PT, annotation_clip=False)
    axis.set_title(title, fontsize=TITLE_PT, pad=4 + line_height_pt * len(stats_lines))


def draw_boxes(axis, groups, labels, colours, point_size=4, show_n=True):
    # Box plots with the individual points jittered on top; the seed fixes the jitter.
    boxes = axis.boxplot(groups, patch_artist=True, widths=0.62, showfliers=False,
                         medianprops=dict(color="black", lw=0.8), whiskerprops=dict(lw=0.5),
                         capprops=dict(lw=0.5), boxprops=dict(lw=0.5))
    for box, colour in zip(boxes["boxes"], colours):
        box.set_facecolor(colour)
        box.set_alpha(0.75)
    random_generator = np.random.default_rng(1)
    for group_index, values in enumerate(groups):
        x_positions = random_generator.normal(group_index + 1, 0.07, len(values))
        axis.scatter(x_positions, values, s=point_size, c="k", alpha=0.55, lw=0, zorder=3)
    tick_labels = []
    for label, values in zip(labels, groups):
        if show_n:
            tick_labels.append(f"{label}\n(n={len(values)})")
        else:
            tick_labels.append(label)
    axis.set_xticks(range(1, len(labels) + 1))
    axis.set_xticklabels(tick_labels, fontsize=SMALL_PT)


def load_states(folder, cell_type):
    # EcoTyper state abundances, returned as samples x states.
    table = pd.read_csv(f"{folder}/{cell_type}/state_abundances.txt", sep="\t", index_col=0)
    table.columns = table.columns.astype(str)
    return table.T


def subtype_colours(subtypes):
    colours = []
    for subtype in subtypes:
        if VERTICAL_GROWTH[subtype]:
            colours.append(VERTICAL_COLOUR)
        else:
            colours.append(RADIAL_COLOUR)
    return colours


# GeoMx annotations shared by Figs. 7 and 8.
geomx_metadata = pd.read_csv(f"{BASE}/data/2025/2025-01-27_GeoMxmetadata.txt", sep="\t")
geomx_metadata["ID"] = geomx_metadata["ID"].astype(str)
geomx_metadata = geomx_metadata.set_index("ID")
geomx_metadata["AOI"] = geomx_metadata["Segment_Tags"].replace({"neg/neg": "neg_neg"})
geomx_assignment = pd.read_csv(f"{RECOVERY}/Ecotypes/ecotype_assignment.txt", sep="\t", index_col=0)
geomx_assignment.index = geomx_assignment.index.astype(str)
geomx_assignment["AOI"] = geomx_assignment["Segment_Tags"].replace({"neg/neg": "neg_neg"})


# ================================================================== Fig. 6a  TCGA-SKCM composition
def fig6a():
    # Step 1. Cell-type fractions per tumour, rescaled to sum to 1, ranked by melanoma fraction.
    fractions = pd.read_csv(f"{BASE}/data/2026-01-28_TCGA_SKCM_fraction.txt", sep="\t", index_col=0)
    fractions = fractions[CELL_TYPE_ORDER].astype(float)
    fractions = fractions.div(fractions.sum(axis=1), axis=0)
    fractions = fractions.sort_values("Melanoma")

    # Step 2. Stacked bars, one per tumour (498 x 140 pt).
    figure = new_figure(498, 140)
    axis = figure.add_subplot()
    positions = np.arange(fractions.shape[0])
    bottom = np.zeros(fractions.shape[0])
    for cell_type in CELL_TYPE_ORDER:
        axis.bar(positions, fractions[cell_type].values, bottom=bottom, width=1.0,
                 color=CELL_COLOUR[cell_type], label=DISPLAY_NAME[cell_type], linewidth=0)
        bottom = bottom + fractions[cell_type].values

    # Step 3. Labels and legend.
    axis.set_xlim(-0.5, fractions.shape[0] - 0.5)
    axis.set_ylim(0, 1)
    axis.set_title("Cell-type composition of TCGA-SKCM tumors", fontsize=TITLE_PT)
    axis.set_xlabel(f"Tumors ranked by melanoma fraction (n={fractions.shape[0]})", fontsize=TITLE_PT)
    axis.set_ylabel("Cell-type fraction", fontsize=TITLE_PT)
    axis.legend(ncol=9, loc="upper center", bbox_to_anchor=(0.5, -0.28), frameon=False,
                handlelength=1, columnspacing=1.2, fontsize=SMALL_PT)
    save_panel(figure, "fig6a_tcga_composition")


# ================================================================== Fig. 6b  melanoma cell states
def fig6b():
    # Step 1. EcoTyper heat-map data for the melanoma compartment and the state of each tumour.
    heatmap_data = pd.read_csv(f"{DISCOVERY}/Melanoma/heatmap_data.txt", sep="\t", index_col=0)
    heatmap_data.columns = heatmap_data.columns.astype(str)
    gene_info = pd.read_csv(f"{DISCOVERY}/Melanoma/gene_info.txt", sep="\t", index_col=0)
    top_annotation = pd.read_csv(f"{DISCOVERY}/Melanoma/heatmap_top_ann.txt", sep="\t", index_col=0)
    top_annotation.index = top_annotation.index.astype(str)
    melanoma_states = load_states(DISCOVERY, "Melanoma")

    # Step 2. Tumours ordered by state (S01 first), then by S01 abundance; genes grouped by state.
    tumors = [t for t in heatmap_data.columns if t in melanoma_states.index]
    s01_tumors = [t for t in tumors if top_annotation.loc[t, "State"] == "S01"]
    s02_tumors = [t for t in tumors if top_annotation.loc[t, "State"] != "S01"]
    s01_tumors = sorted(s01_tumors, key=lambda t: -melanoma_states.loc[t, "S01"])
    s02_tumors = sorted(s02_tumors, key=lambda t: -melanoma_states.loc[t, "S01"])
    tumor_order = s01_tumors + s02_tumors
    s01_genes = [g for g in heatmap_data.index if gene_info.loc[g, "State"] == "S01"]
    s02_genes = [g for g in heatmap_data.index if gene_info.loc[g, "State"] == "S02"]
    gene_order = s01_genes + s02_genes
    values = np.clip(heatmap_data.loc[gene_order, tumor_order].values, -2.5, 2.5)   # clip extremes for colour range
    n_tumors = len(tumor_order)

    # Step 3. Layout (255 x 220 pt): three annotation bars over the heat map, gene labels on the right.
    figure = new_figure(255, 220)
    grid = GridSpec(4, 2, figure=figure, height_ratios=[0.5, 0.5, 0.5, 12], width_ratios=[30, 6.5], hspace=0.04, wspace=0.0)

    def annotation_bar(row, bar_values, colour_map, label):
        bar_axis = figure.add_subplot(grid[row, 0])
        vector_heatmap(bar_axis, bar_values[None, :], colour_map)
        bar_axis.set_xticks([])
        bar_axis.set_yticks([])
        for spine in bar_axis.spines.values():
            spine.set_visible(False)
        bar_axis.text(-0.01, 0.5, label, transform=bar_axis.transAxes, ha="right", va="center", fontsize=SMALL_PT)
        return bar_axis

    state_codes = np.array([0] * len(s01_tumors) + [1] * len(s02_tumors))
    first_bar = annotation_bar(0, state_codes, ListedColormap(["#E8336D", "#4BA3A3"]), "State")
    first_bar.set_title("Melanoma cell states (TCGA-SKCM)", fontsize=TITLE_PT, pad=3)
    sample_type_codes = []
    for tumor in tumor_order:
        sample_type = str(top_annotation.loc[tumor, "SAMPLE_TYPE"])
        if sample_type == "Primary":
            sample_type_codes.append(0)
        elif sample_type == "Metastasis":
            sample_type_codes.append(1)
        else:
            sample_type_codes.append(2)
    annotation_bar(1, np.array(sample_type_codes), ListedColormap([PRIMARY_COLOUR, METASTASIS_COLOUR, MISSING_COLOUR]), "Sample type")
    annotation_bar(2, melanoma_states.loc[tumor_order, "S01"].values, "viridis", "S01 abundance")

    # Step 4. Heat map with state blocks separated by black lines.
    heatmap_axis = figure.add_subplot(grid[3, 0])
    # The 951 x 473 expression heat map is the one raster element (author decision, 30 Oct 2026): as vector
    # rectangles it made Fig. 6 about 15 MB. imshow embeds it at its native resolution (one image pixel per
    # gene x tumor cell, no smoothing), so it stays sharp at any zoom; all other Fig. 6b elements stay vector.
    image = heatmap_axis.imshow(values, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5, interpolation="none")   # "none": no resampling in PDF
    heatmap_axis.axhline(len(s01_genes) - 0.5, color="k", lw=0.8)
    heatmap_axis.axvline(len(s01_tumors) - 0.5, color="k", lw=0.8)
    heatmap_axis.set_xticks([])
    heatmap_axis.set_yticks([])
    heatmap_axis.set_xlabel(f"Tumors (n={n_tumors}), ordered by state and S01 abundance", fontsize=TITLE_PT)
    heatmap_axis.text(-0.02, 1 - len(s01_genes) / 2 / len(gene_order), f"S01\npigmented\n({len(s01_genes)} genes)",
                      transform=heatmap_axis.transAxes, ha="right", va="center", fontsize=SMALL_PT, color="#C2185B")
    heatmap_axis.text(-0.02, (len(s02_genes) / 2) / len(gene_order), f"S02\nde-diff.\n({len(s02_genes)} genes)",
                      transform=heatmap_axis.transAxes, ha="right", va="center", fontsize=SMALL_PT, color="#2E7D7D")

    # Step 5. Marker-gene labels on the right, joined to their rows by leader lines; labels are spaced
    # at least min_gap rows apart so they never overlap.
    s01_markers = ["OCA2", "PMEL", "GPR143", "TRPM1", "TSPAN10", "MFSD12", "MLANA", "TYR", "MITF"]
    s02_markers = ["AXL", "NGFR", "SOX9", "WNT5A", "TGFBI", "SERPINE1", "IL6ST", "AHR", "LEF1", "CD68", "AEBP1", "SERPINE2"]
    label_axis = figure.add_subplot(grid[3, 1])
    label_axis.set_xlim(0, 1)
    label_axis.set_ylim(len(gene_order) - 0.5, -0.5)
    label_axis.axis("off")
    labelled_genes = []
    for gene in s01_markers + s02_markers:
        if gene in gene_order:
            labelled_genes.append((gene_order.index(gene), gene))
    labelled_genes = sorted(labelled_genes)
    label_rows = np.array([row for row, _ in labelled_genes], dtype=float)
    min_gap = len(gene_order) / 18          # about 8 pt between labels at this panel height
    for index in range(1, len(label_rows)):
        if label_rows[index] - label_rows[index - 1] < min_gap:
            label_rows[index] = label_rows[index - 1] + min_gap
    # Each leader line is an explicit elbow in heat-map axes coordinates: out of the gene's own row,
    # diagonally to the label's height, then into the label's left edge, so every line ends at its label.
    n_genes = len(gene_order)
    for (gene_row, gene), label_row in zip(labelled_genes, label_rows):
        row_y = 1 - (gene_row + 0.5) / n_genes
        label_y = 1 - (label_row + 0.5) / n_genes
        line_x = [1.0, 1.02, 1.08, 1.10]
        line_y = [row_y, row_y, label_y, label_y]
        heatmap_axis.plot(line_x, line_y, transform=heatmap_axis.transAxes, color="grey", lw=0.4, clip_on=False)
        heatmap_axis.text(1.11, label_y, italic(gene), transform=heatmap_axis.transAxes,
                          fontsize=SMALL_PT, va="center", ha="left")

    # Step 6. Colour bar and legend for the annotation bars.
    colorbar = keep_colorbar_vector(figure.colorbar(image, ax=heatmap_axis, location="bottom", fraction=0.03, pad=0.03, shrink=0.3, anchor=(1.0, 0)))
    colorbar.set_label("Expression (z)", fontsize=TITLE_PT)
    colorbar.ax.tick_params(labelsize=SMALL_PT)
    legend_handles = [Patch(color="#E8336D", label="S01"), Patch(color="#4BA3A3", label="S02"),
                      Patch(color=PRIMARY_COLOUR, label="Primary"), Patch(color=METASTASIS_COLOUR, label="Metastasis"),
                      Patch(color=MISSING_COLOUR, label="Not annotated")]
    figure.legend(handles=legend_handles, loc="outside lower left", ncol=5, frameon=False, fontsize=SMALL_PT,
                  handlelength=1, columnspacing=0.8, handletextpad=0.4)
    save_panel(figure, "fig6b_tcga_melanoma_states")


# ================================================================== Fig. 6c  S01 versus signatures
PIGMENTATION_GENES = ["PMEL", "MLANA", "TYR", "TYRP1", "DCT", "MITF", "OCA2", "GPR143", "TRPM1", "SLC45A2", "MFSD12", "TSPAN10"]
IFN_GENES = ["STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11", "GBP1", "IDO1", "HLA-DRA", "CD8A", "GZMB", "PRF1", "CD274", "TAP1"]


def fig6c():
    # Step 1. Read only the signature genes from the large TCGA TPM file.
    wanted_genes = set(PIGMENTATION_GENES + IFN_GENES)
    gene_rows = []
    with open(f"{BASE}/data/2025/2025-01-29_SKCM_coding_symbol_TPM.txt") as tpm_file:
        header = tpm_file.readline().rstrip("\n").split("\t")
        for line in tpm_file:
            gene = line.split("\t", 1)[0]
            if gene in wanted_genes:
                gene_rows.append(line.rstrip("\n").split("\t"))
    tpm = pd.DataFrame([row[1:] for row in gene_rows], index=[row[0] for row in gene_rows], columns=header[1:]).astype(float)

    # Step 2. Signature scores = mean z-score of log2(TPM + 1) across the signature genes.
    melanoma_states = load_states(DISCOVERY, "Melanoma")
    tumors = [t for t in melanoma_states.index if t in tpm.columns]
    log_tpm = np.log2(tpm[tumors] + 1)
    z_scores = log_tpm.sub(log_tpm.mean(axis=1), axis=0).div(log_tpm.std(axis=1) + 1e-9, axis=0)
    pigmentation_score = z_scores.loc[[g for g in PIGMENTATION_GENES if g in z_scores.index]].mean(axis=0)
    ifn_score = z_scores.loc[[g for g in IFN_GENES if g in z_scores.index]].mean(axis=0)
    s01_abundance = melanoma_states.loc[tumors, "S01"]

    # Step 3. Two stacked scatter plots (219 x 220 pt).
    figure = new_figure(219, 220)
    axes = figure.subplots(2, 1)
    panels = [(axes[0], pigmentation_score, "Pigmentation/MITF\nsignature (z)"),
              (axes[1], ifn_score, "IFN-γ/inflammation\nsignature (z)")]
    for axis, score, y_label in panels:
        rho, p_value = stats.spearmanr(s01_abundance, score)
        print(f"6c {y_label.splitlines()[0]}: rho={rho:.2f} P={p_value:.2g}")
        axis.scatter(s01_abundance, score, c=s01_abundance, cmap="RdBu", vmin=0.15, vmax=0.85, s=5, edgecolor="k", lw=0.15)
        slope, intercept = np.polyfit(s01_abundance, score, 1)
        line_x = np.linspace(0, 1, 5)
        axis.plot(line_x, slope * line_x + intercept, color="grey", ls="--", lw=0.7)
        axis.set_ylabel(y_label, fontsize=TITLE_PT)
        axis.set_xlim(-0.02, 1.02)
        stats_text = f"Spearman ρ = {rho:.2f}\n$\\it{{P}}$ = {format_p(p_value)}"
        white_backing = dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1)
        if rho > 0:
            axis.text(0.03, 0.97, stats_text, transform=axis.transAxes, fontsize=SMALL_PT, va="top", ha="left", bbox=white_backing)
        else:
            axis.text(0.97, 0.97, stats_text, transform=axis.transAxes, fontsize=SMALL_PT, va="top", ha="right", bbox=white_backing)
    axes[0].set_title("S01 abundance versus expression signatures", fontsize=TITLE_PT)
    axes[0].tick_params(labelbottom=False)
    axes[1].set_xlabel(f"Melanoma-S01 abundance (n={len(tumors)} tumors)", fontsize=TITLE_PT)
    save_panel(figure, "fig6c_s01_vs_signatures")


# ================================================================== Figs. 6d and 7a  EcoTyper ecotype heat maps
def ecotype_heatmap(run_folder, assignment, annotation_tracks, title, name, width_pt, height_pt):
    """Heat map of cell-state abundances for the samples assigned to an ecotype.

    Rows: the cell states that define E1 and E2 (from the discovery run's ecotypes.txt), E1 first.
    Columns: assigned samples, E1 first; within an ecotype, ordered by that ecotype's abundance.
    Values: raw EcoTyper state abundance (no rescaling), viridis colour scale from 0 to the maximum.
    annotation_tracks: list of (track label, values per sample, {category: colour}).
    """
    # Step 1. Ecotype definitions and the state abundance of each defining state in each sample.
    ecotype_states = pd.read_csv(f"{DISCOVERY}/Ecotypes/ecotypes.txt", sep="\t")
    ecotype_states = ecotype_states.sort_values(["Ecotype"], kind="stable")
    ecotype_abundance = pd.read_csv(f"{run_folder}/Ecotypes/ecotype_abundance.txt", sep="\t", index_col=0)
    ecotype_abundance.columns = ecotype_abundance.columns.astype(str)

    # Step 2. Column order: E1 samples then E2 samples, each by decreasing abundance of its own ecotype.
    sample_order = []
    for ecotype in ["E1", "E2"]:
        samples = assignment.index[assignment["Ecotype"] == ecotype].tolist()
        samples = sorted(samples, key=lambda s: -ecotype_abundance.loc[ecotype, s])
        sample_order = sample_order + samples
    n_e1_samples = int((assignment["Ecotype"] == "E1").sum())

    # Step 3. Matrix of state abundances (rows = ecotype states, columns = samples).
    row_values = []
    row_cell_types = []
    row_states = []
    row_ecotypes = []
    state_tables = {}
    for _, state_row in ecotype_states.iterrows():
        cell_type = state_row["CellType"]
        if cell_type not in state_tables:
            state_tables[cell_type] = load_states(run_folder, cell_type)
        row_values.append(state_tables[cell_type].loc[sample_order, state_row["State"]].values)
        row_cell_types.append(cell_type.replace(".", "-"))      # EcoTyper writes T.CD8; the palette uses T-CD8
        row_states.append(state_row["State"])
        row_ecotypes.append(state_row["Ecotype"])
    matrix = np.array(row_values, dtype=float)
    n_e1_rows = row_ecotypes.count("E1")
    print(f"{name}: {matrix.shape[0]} states x {matrix.shape[1]} samples (E1 {n_e1_samples}), max abundance {matrix.max():.2f}")

    # Step 4. Layout: annotation tracks, the heat map, then a legend row. The legend row is split into one
    # column per annotation track (each sized from its longest text, so blocks cannot overlap) and a second
    # line for the colour bar. Legends are part of the layout, so the heat map shrinks to make room for them.
    n_tracks = len(annotation_tracks)
    figure = new_figure(width_pt, height_pt)
    grid = GridSpec(n_tracks + 2, 1, figure=figure, height_ratios=[0.55] * n_tracks + [11, 10], hspace=0.05)

    # Step 5. Annotation tracks above the heat map, each labelled on the left.
    for track_index, (track_label, track_values, track_colours) in enumerate(annotation_tracks):
        track_axis = figure.add_subplot(grid[track_index, 0])
        categories = list(track_colours.keys())
        codes = np.array([categories.index(track_values.loc[s]) for s in sample_order])
        vector_heatmap(track_axis, codes[None, :], ListedColormap([track_colours[c] for c in categories]),
                       vmin=0, vmax=len(categories) - 1)
        track_axis.axvline(n_e1_samples - 0.5, color="white", lw=0.8)
        track_axis.set_xticks([])
        track_axis.set_yticks([])
        for spine in track_axis.spines.values():
            spine.set_visible(False)
        track_axis.text(-0.01, 0.5, track_label, transform=track_axis.transAxes, ha="right", va="center", fontsize=SMALL_PT)
        if track_index == 0:
            track_axis.set_title(title, fontsize=TITLE_PT, pad=3)

    # Step 6. Heat map. Rows are identified by the colour strips on the left (Step 7), not by text labels.
    heatmap_axis = figure.add_subplot(grid[n_tracks, 0])
    image = vector_heatmap(heatmap_axis, matrix, "viridis", vmin=0, vmax=matrix.max())
    heatmap_axis.axhline(n_e1_rows - 0.5, color="white", lw=1.0)
    heatmap_axis.axvline(n_e1_samples - 0.5, color="white", lw=1.0)
    heatmap_axis.set_xticks([])
    heatmap_axis.set_yticks([])
    heatmap_axis.set_xlabel(f"Assigned samples (n={len(sample_order)})", fontsize=TITLE_PT)

    # Step 7. Three row strips attached to the left edge of the heat map, as in the EcoTyper output:
    # Ecotype, Cell type and State (left to right). Each strip is labelled underneath, rotated 90 degrees.
    state_colours = {"S01": "#E8336D", "S02": "#4BA3A3", "S03": "#B07AA1"}   # S01/S02 match Fig. 6b
    strip_specs = [
        ("Ecotype", row_ecotypes, {"E1": E1_COLOUR, "E2": E2_COLOUR}),
        ("Cell type", row_cell_types, {c: CELL_COLOUR[c] for c in CELL_TYPE_ORDER if c in row_cell_types}),
        ("State", row_states, {s: state_colours[s] for s in sorted(set(row_states))}),
    ]
    strip_width = 0.022          # fraction of the heat-map width per strip
    strip_gap = 0.006
    for strip_index, (strip_label, strip_values, strip_colours) in enumerate(strip_specs):
        categories = list(strip_colours.keys())
        left = -(len(strip_specs) - strip_index) * (strip_width + strip_gap)
        strip_axis = heatmap_axis.inset_axes([left, 0, strip_width, 1])
        codes = np.array([categories.index(v) for v in strip_values])[:, None]
        vector_heatmap(strip_axis, codes, ListedColormap([strip_colours[c] for c in categories]),
                       vmin=0, vmax=len(categories) - 1)
        strip_axis.axhline(n_e1_rows - 0.5, color="white", lw=1.0)
        strip_axis.set_xticks([])
        strip_axis.set_yticks([])
        for spine in strip_axis.spines.values():
            spine.set_visible(False)
        strip_axis.text(0.5, -0.03, strip_label, transform=strip_axis.transAxes, rotation=90,
                        ha="center", va="top", fontsize=SMALL_PT)

    # Step 8. Legend area, two lines. Line 1: one block per annotation track. Line 2: Cell type, State and
    # the abundance colour bar. Column widths come from the longest label or title in each block (about
    # 0.55 x font size per character, plus the colour key); long lists use three columns.
    def block_width(title_text, categories, n_columns):
        longest_label_pt = max(len(str(c)) for c in categories) * SMALL_PT * 0.55 + 10
        title_pt = len(title_text) * TITLE_PT * 0.55
        return max(title_pt, longest_label_pt * n_columns)

    line1_specs = []
    for track_label, _, track_colours in annotation_tracks:
        n_columns = 1 if len(track_colours) <= 4 else 3
        line1_specs.append((track_label, track_colours, n_columns))
    cell_type_names = {c: DISPLAY_NAME[c] for c in strip_specs[1][2]}
    line2_specs = [("Cell type", {cell_type_names[c]: colour for c, colour in strip_specs[1][2].items()}, 3),
                   ("State", strip_specs[2][2], 1)]

    legend_area = grid[n_tracks + 1, 0].subgridspec(2, 1, height_ratios=[1.15, 1], hspace=0.25)

    def draw_legend_line(subplot_spec, specs, extra_width_pt):
        widths = [block_width(title_text, categories, n_columns) for title_text, categories, n_columns in specs]
        line_grid = subplot_spec.subgridspec(1, len(specs) + 1, width_ratios=widths + [extra_width_pt], wspace=0.15)
        for column, (title_text, categories, n_columns) in enumerate(specs):
            legend_axis = figure.add_subplot(line_grid[0, column])
            legend_axis.axis("off")
            handles = [Patch(color=colour, label=category) for category, colour in categories.items()]
            block = legend_axis.legend(handles=handles, title=title_text, loc="upper left", bbox_to_anchor=(0, 1),
                                       frameon=False, fontsize=SMALL_PT, title_fontsize=TITLE_PT, ncol=n_columns,
                                       handlelength=0.9, handletextpad=0.4, columnspacing=0.6, borderpad=0, labelspacing=0.25)
            block._legend_box.align = "left"
        return line_grid

    draw_legend_line(legend_area[0, 0], line1_specs, 20)
    line2_grid = draw_legend_line(legend_area[1, 0], line2_specs, 70)

    # Colour bar in the last column of line 2, 60 pt wide.
    colorbar_holder = figure.add_subplot(line2_grid[0, len(line2_specs)])
    colorbar_holder.axis("off")
    colorbar_axis = colorbar_holder.inset_axes([0.05, 0.55, 0.85, 0.12])
    colorbar = keep_colorbar_vector(figure.colorbar(image, cax=colorbar_axis, orientation="horizontal"))
    colorbar.set_label("State abundance", fontsize=TITLE_PT, labelpad=1)
    colorbar.ax.xaxis.set_label_position("top")
    colorbar.ax.tick_params(labelsize=SMALL_PT, length=1.5, pad=1)
    save_panel(figure, name)


def fig6d():
    # TCGA-SKCM discovery: 358 assigned tumours with sample type, aneuploidy and TMB annotation.
    assignment = pd.read_csv(f"{DISCOVERY}/Ecotypes/ecotype_assignment.txt", sep="\t", index_col=0)
    assignment.index = assignment.index.astype(str)
    sample_type = assignment["SAMPLE_TYPE"].fillna("Not annotated").replace({"nan": "Not annotated"})
    aneuploidy_names = {"ANEUPLOIDY<7": "<7", "<7ANEUPLOIDY<20": "7–20", "ANEUPLOIDY>20": ">20"}
    aneuploidy = assignment["ANEUPLOIDY_Category"].map(aneuploidy_names).fillna("Not annotated")
    tmb_names = {"TMB<10": "<10", "TMB>=10": "≥10"}
    tmb = assignment["TMB_Category"].map(tmb_names).fillna("Not annotated")
    ecotype = assignment["Ecotype"]
    tracks = [
        ("Ecotype", ecotype, {"E1": E1_COLOUR, "E2": E2_COLOUR}),
        ("Sample type", sample_type, {"Primary": PRIMARY_COLOUR, "Metastasis": METASTASIS_COLOUR, "Not annotated": MISSING_COLOUR}),
        ("Aneuploidy", aneuploidy, {"<7": "#C7E9C0", "7–20": "#74C476", ">20": "#238B45", "Not annotated": MISSING_COLOUR}),
        ("TMB (mut/Mb)", tmb, {"<10": "#C6DBEF", "≥10": "#2171B5", "Not annotated": MISSING_COLOUR}),
    ]
    ecotype_heatmap(DISCOVERY, assignment, tracks, "TCGA-SKCM ecotypes", "fig6d_tcga_ecotype_heatmap", 262, 270)


def fig7a():
    # GeoMx recovery: 167 assigned AOIs with segment, tissue and subtype annotation.
    assignment = geomx_assignment.copy()
    segment = assignment["Segment_Tags"].map({"S100": "S100+", "CD45": "CD45+", "neg/neg": "S100−/CD45−"})
    tissue = assignment["Pathophysiology"]
    subtype = assignment["Path_sub_type"]
    tracks = [
        ("Ecotype", assignment["Ecotype"], {"E1": E1_COLOUR, "E2": E2_COLOUR}),
        ("Segment", segment, {"S100+": "#619CFF", "CD45+": "#F8766D", "S100−/CD45−": "#00BA38"}),
        ("Tissue", tissue, {"Nevi": "#BDBDBD", "Melanoma": "#525252"}),
        ("Subtype", subtype, {"Junc": "#E5F5E0", "Comp": "#A1D99B", "ID": "#41AB5D",
                              "LMM": "#C6DBEF", "SITU": "#6BAED6", "SSM": "#2171B5",
                              "NOD": "#FCBBA1", "AMEL": "#FB6A4A", "MET": "#A50F15"}),
    ]
    ecotype_heatmap(RECOVERY, assignment, tracks, "GeoMx ecotypes (recovered)", "fig7a_geomx_ecotype_heatmap", 320, 255)


# ================================================================== Fig. 6e  ecotype associations
def fig6e():
    # Step 1. TCGA ecotype assignments with sample type, aneuploidy and TMB.
    assignment = pd.read_csv(f"{DISCOVERY}/Ecotypes/ecotype_assignment.txt", sep="\t", index_col=0)
    assignment["ANEUPLOIDY_SCORE"] = pd.to_numeric(assignment["ANEUPLOIDY_SCORE"], errors="coerce")
    assignment["TMB_NONSYNONYMOUS"] = pd.to_numeric(assignment["TMB_NONSYNONYMOUS"], errors="coerce")

    # Step 2. Layout (212 x 235 pt): sample-type bars and two box plots.
    figure = new_figure(212, 235)
    axes = figure.subplots(1, 3, width_ratios=[1.25, 1, 1])

    # Step 3. Sample type by ecotype; only tumours with a sample-type annotation enter Fisher's test.
    sample_type_axis = axes[0]
    annotated = assignment[assignment["SAMPLE_TYPE"].isin(["Primary", "Metastasis"])]
    counts = pd.crosstab(annotated["Ecotype"], annotated["SAMPLE_TYPE"])[["Primary", "Metastasis"]]
    odds_ratio, fisher_p = stats.fisher_exact(counts.values)
    print(f"6e Fisher: n={len(annotated)} OR(met odds, E1 vs E2)={1 / odds_ratio:.2f} P={fisher_p:.2g}")
    percent = counts.div(counts.sum(axis=1), axis=0) * 100
    sample_type_axis.bar([0, 1], percent["Metastasis"], color=METASTASIS_COLOUR, width=0.65, label="Metastasis")
    sample_type_axis.bar([0, 1], percent["Primary"], bottom=percent["Metastasis"], color=PRIMARY_COLOUR, width=0.65, label="Primary")
    sample_type_axis.set_xticks([0, 1])
    sample_type_axis.set_xticklabels([f"E1\nn={counts.loc['E1'].sum()}", f"E2\nn={counts.loc['E2'].sum()}"], fontsize=SMALL_PT)
    sample_type_axis.set_ylabel("% of tumors", fontsize=TITLE_PT)
    sample_type_axis.set_ylim(0, 100)
    title_with_stats(sample_type_axis, "Sample type", [f"Fisher $\\it{{P}}$ = {format_p(fisher_p)}"])
    sample_type_axis.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), frameon=False, ncol=1, fontsize=SMALL_PT, handlelength=0.9)

    # Step 4. Aneuploidy and TMB by ecotype (Mann-Whitney).
    box_panels = [(axes[1], "ANEUPLOIDY_SCORE", "Aneuploidy", "Aneuploidy score", False),
                  (axes[2], "TMB_NONSYNONYMOUS", "TMB", "TMB (mut/Mb)", True)]
    for axis, column, title, y_label, log_scale in box_panels:
        e1_values = assignment[assignment["Ecotype"] == "E1"][column].dropna()
        e2_values = assignment[assignment["Ecotype"] == "E2"][column].dropna()
        mann_whitney_p = stats.mannwhitneyu(e1_values, e2_values).pvalue
        print(f"6e {title}: P={mann_whitney_p:.2g}")
        draw_boxes(axis, [e1_values.values, e2_values.values], ["E1", "E2"], [E1_COLOUR, E2_COLOUR], point_size=2, show_n=False)
        if log_scale:
            axis.set_yscale("log")
        axis.set_ylabel(y_label, fontsize=TITLE_PT)
        title_with_stats(axis, title, [f"$\\it{{P}}$ = {format_p(mann_whitney_p)} (n.s.)"])
    save_panel(figure, "fig6e_tcga_associations")


# ================================================================== Fig. 7b  ecotype by segment
def fig7b():
    # Step 1. Percentage of melanoma AOIs assigned to E1 or E2 in each segment.
    melanoma_aois = geomx_assignment[geomx_assignment["Pathophysiology"] == "Melanoma"]
    segments = [("S100", "Tumor\nS100+"), ("CD45", "Immune\nCD45+"), ("neg_neg", "Stromal\nS100−\nCD45−")]
    segment_keys = [key for key, _ in segments]
    counts = pd.crosstab(melanoma_aois["AOI"], melanoma_aois["Ecotype"]).reindex(segment_keys)
    percent = counts.div(counts.sum(axis=1), axis=0) * 100
    print("7b E2 % by segment:", percent["E2"].round(1).to_dict())

    # Step 2. Stacked bars (154 x 225 pt).
    figure = new_figure(154, 225)
    axis = figure.add_subplot()
    positions = np.arange(len(segments))
    axis.bar(positions, percent["E2"], color=E2_COLOUR, width=0.65, label="E2", edgecolor="k", lw=0.4)
    axis.bar(positions, percent["E1"], bottom=percent["E2"], color=E1_COLOUR, width=0.65, label="E1", edgecolor="k", lw=0.4)
    for position, segment in zip(positions, segment_keys):
        axis.text(position, percent.loc[segment, "E2"] / 2, f"{percent.loc[segment, 'E2']:.1f}%", ha="center", va="center", fontsize=SMALL_PT)
        axis.text(position, 101, f"n={int(counts.loc[segment].sum())}", ha="center", va="bottom", fontsize=SMALL_PT)

    # Step 3. Labels and legend.
    axis.set_xticks(positions)
    axis.set_xticklabels([label for _, label in segments], fontsize=SMALL_PT)
    axis.set_ylim(0, 108)
    axis.set_yticks([0, 25, 50, 75, 100])
    axis.set_ylabel("% of melanoma AOIs", fontsize=TITLE_PT)
    axis.set_title("Ecotype by segment", fontsize=TITLE_PT)
    axis.legend(loc="upper center", bbox_to_anchor=(0.5, -0.24), ncol=2, frameon=False, fontsize=SMALL_PT, handlelength=1)
    save_panel(figure, "fig7b_ecotype_by_segment")


# ================================================================== Fig. 7c  niche coupling (62 melanoma regions)
def fig7c():
    # Niche coupling in melanoma regions (author decision, 30 Oct 2026): every melanoma AOI is labelled with
    # its higher-abundance ecotype, including AOIs below EcoTyper's assignment threshold. For the 167 AOIs
    # EcoTyper did assign, this label always equals EcoTyper's call; the threshold only removes low-confidence
    # AOIs. This lets the test use all 62 melanoma regions with three profiled segments instead of 20.

    # Step 1. Higher-abundance ecotype for every melanoma AOI in the recovery input.
    ecotype_abundance = pd.read_csv(f"{RECOVERY}/Ecotypes/ecotype_abundance.txt", sep="\t", index_col=0).T
    ecotype_abundance.index = ecotype_abundance.index.astype(str)
    labels = ecotype_abundance.idxmax(axis=1).rename("ecotype").to_frame()
    labels = labels.join(geomx_metadata[["Scan_ID", "ROI..Label.", "Segment_Tags", "Pathophysiology"]])
    labels = labels[labels["Pathophysiology"] == "Melanoma"]
    labels["region"] = labels["Scan_ID"] + "|" + labels["ROI..Label."]

    # Step 2. One row per region, one column per segment.
    compartments = ["S100", "CD45", "neg/neg"]
    by_region = labels.pivot_table(index="region", columns="Segment_Tags", values="ecotype", aggfunc="first")
    complete_regions = by_region.dropna(subset=compartments)[compartments]
    n_regions = len(complete_regions)

    # Step 3. All-three concordance versus the rate expected if segments were labelled independently,
    # using the overall E1 frequency among melanoma AOIs.
    e1_fraction = (labels["ecotype"] == "E1").mean()
    expected_concordance = e1_fraction ** 3 + (1 - e1_fraction) ** 3
    is_concordant = complete_regions.nunique(axis=1) == 1
    n_concordant = int(is_concordant.sum())
    binomial_p = stats.binomtest(n_concordant, n_regions, expected_concordance, alternative="greater").pvalue
    print(f"7c: {n_concordant}/{n_regions} concordant ({100 * n_concordant / n_regions:.1f}%), "
          f"expected {100 * expected_concordance:.1f}%, binomial P={binomial_p:.2g}")

    # Step 4. Pairwise agreement, the agreement expected by chance for that pair, and Cohen's kappa.
    pair_specs = [("S100", "CD45", "S100+ vs\nCD45+"), ("S100", "neg/neg", "S100+ vs\nstromal"), ("CD45", "neg/neg", "CD45+ vs\nstromal")]
    pair_results = []
    for first, second, pair_label in pair_specs:
        pair_table = by_region[[first, second]].dropna()
        first_labels = pair_table[first].values
        second_labels = pair_table[second].values
        observed = np.mean(first_labels == second_labels)
        chance = 0.0
        for ecotype in ["E1", "E2"]:
            chance = chance + np.mean(first_labels == ecotype) * np.mean(second_labels == ecotype)
        kappa = (observed - chance) / (1 - chance)
        pair_results.append((pair_label, len(pair_table), 100 * observed, 100 * chance, kappa))
        print(f"7c {pair_label.replace(chr(10), ' ')}: n={len(pair_table)}, agreement {100 * observed:.1f}%, chance {100 * chance:.1f}%, kappa {kappa:.2f}")

    # Step 5. Layout (230 x 150 pt): label matrix on the left, pairwise agreement on the right.
    figure = new_figure(230, 150)
    matrix_axis, agreement_axis = figure.subplots(1, 2, width_ratios=[1.1, 1])

    # Step 6. Label matrix: regions sorted so identical patterns sit together; thin white lines
    # separate the individual regions.
    is_e1 = (complete_regions == "E1").astype(int).values
    region_order = np.lexsort((is_e1[:, 2], is_e1[:, 1], is_e1[:, 0]))[::-1]
    vector_heatmap(matrix_axis, is_e1[region_order], ListedColormap([E2_COLOUR, E1_COLOUR]), vmin=0, vmax=1)
    for boundary in range(1, n_regions):
        matrix_axis.axhline(boundary - 0.5, color="white", lw=0.15)
    matrix_axis.set_xticks(range(3))
    matrix_axis.set_xticklabels(["S100+", "CD45+", "S100−/\nCD45−"], fontsize=SMALL_PT)
    matrix_axis.set_yticks([])
    matrix_axis.set_ylabel(f"Melanoma regions (n={n_regions})", fontsize=TITLE_PT)
    for spine in matrix_axis.spines.values():
        spine.set_visible(True)
    title_with_stats(matrix_axis, "Ecotype per compartment",
                     [f"{n_concordant}/{n_regions} concordant ({100 * n_concordant / n_regions:.0f}%;",
                      f"expected {100 * expected_concordance:.0f}%)"])
    matrix_axis.legend(handles=[Patch(color=E1_COLOUR, label="E1"), Patch(color=E2_COLOUR, label="E2")],
                       loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False, fontsize=SMALL_PT, handlelength=0.9)

    # Step 7. Pairwise agreement bars; a dashed tick on each bar marks that pair's chance agreement.
    for index, (pair_label, n_pair, observed_pct, chance_pct, kappa) in enumerate(pair_results):
        agreement_axis.bar(index, observed_pct, color="#8E7CC3", width=0.6)
        agreement_axis.plot([index - 0.36, index + 0.36], [chance_pct, chance_pct], color="black", ls="--", lw=0.6)
        kappa_text = f"κ = {kappa:.2f}".replace("-", "−")
        agreement_axis.text(index, max(observed_pct, chance_pct) + 3, kappa_text, ha="center", fontsize=SMALL_PT)
    agreement_axis.set_xticks(range(len(pair_results)))
    agreement_axis.set_xticklabels([f"{label}\n(n={n})" for label, n, _, _, _ in pair_results], fontsize=SMALL_PT)
    agreement_axis.set_ylim(0, 100)
    agreement_axis.set_ylabel("Pairwise agreement (%)", fontsize=TITLE_PT)
    title_with_stats(agreement_axis, "Pairwise agreement", ["dashed line: chance", ""])
    save_panel(figure, "fig7c_niche_coupling")


# ================================================================== Figs. 7d-f  cell states by subtype
def subtype_groups(states, aoi_type, cell_type, state):
    # State abundance per melanoma subtype within one AOI type, with the Kruskal-Wallis test across
    # subtypes and the Mann-Whitney test of radial versus vertical growth (same tests as the text).
    # AOI set (author decision, 10 Nov 2026): every melanoma AOI of this type, whether or not EcoTyper
    # assigned it to an ecotype. EcoTyper estimates state abundances for all AOIs, so the ecotype
    # threshold is not needed here; until 10 Nov 2026 only the assigned AOIs were used.
    melanoma_aois = geomx_metadata[geomx_metadata["Pathophysiology"] == "Melanoma"]
    aoi_ids = melanoma_aois.index[melanoma_aois["AOI"] == aoi_type]
    table = pd.DataFrame({"value": states[cell_type][state].reindex(aoi_ids),
                          "subtype": melanoma_aois.loc[aoi_ids, "Path_sub_type"]}).dropna()
    subtypes_present = []
    for subtype in SUBTYPES:
        if (table["subtype"] == subtype).sum() >= 3:      # subtypes with fewer than 3 AOIs are not shown
            subtypes_present.append(subtype)
    groups = []
    for subtype in subtypes_present:
        groups.append(table.loc[table["subtype"] == subtype, "value"].values)
    kruskal_p = stats.kruskal(*groups).pvalue
    vertical = table["subtype"].map(VERTICAL_GROWTH)
    radial_vs_vertical_p = stats.mannwhitneyu(table.loc[~vertical, "value"], table.loc[vertical, "value"]).pvalue
    return subtypes_present, groups, kruskal_p, radial_vs_vertical_p


def fig7def():
    # Step 1. Recovered state abundances for every cell type.
    cell_types = ["Melanoma", "Endothelial", "CAF", "Mural_cell", "B_cell", "T.CD8", "NK"]
    states = {}
    for cell_type in cell_types:
        states[cell_type] = load_states(RECOVERY, cell_type)

    # Step 2. Fig. 7d: pigmented Melanoma-S01 in S100+ tumor AOIs (244 x 150 pt).
    subtypes, groups, kruskal_p, _ = subtype_groups(states, "S100", "Melanoma", "S01")
    print("7d medians:", {s: round(float(np.median(g)), 2) for s, g in zip(subtypes, groups)}, "KW P", f"{kruskal_p:.2g}")
    figure = new_figure(244, 150)
    axis = figure.add_subplot()
    draw_boxes(axis, groups, subtypes, subtype_colours(subtypes))
    axis.axhline(0.5, color="grey", ls="--", lw=0.6)
    axis.set_ylabel("Melanoma-S01 abundance", fontsize=TITLE_PT)
    axis.set_xlabel("Histologic subtype (S100+ tumor AOIs)", fontsize=TITLE_PT)
    title_with_stats(axis, "Pigmented melanoma state", [f"Kruskal–Wallis $\\it{{P}}$ = {format_p(kruskal_p)}"])
    axis.legend(handles=[Patch(color=RADIAL_COLOUR, label="Radial growth"), Patch(color=VERTICAL_COLOUR, label="Vertical growth / metastasis")],
                loc="lower left", frameon=False, fontsize=SMALL_PT, handlelength=1)
    save_panel(figure, "fig7d_S01_by_subtype")

    # Step 3. Figs. 7e and 7f: three states each, one row per figure panel (498 x 140 pt).
    panel_specs = [
        ("fig7e_stromal_states", "neg_neg", "Stromal (S100−/CD45−) AOIs",
         [("Endothelial", "S01", "Endothelial-S01 (quiescent vessel)"),
          ("CAF", "S01", "CAF-S01 (myofibroblast/desmoplastic)"),
          ("Mural_cell", "S01", "Mural-S01 (mature pericyte)")]),
        ("fig7f_immune_states", "CD45", "Immune (CD45+) AOIs",
         [("B_cell", "S02", "B cell-S02"),
          ("T.CD8", "S01", "CD8 T-S01 (memory)"),
          ("NK", "S01", "NK-S01")]),
    ]
    for file_name, aoi_type, row_title, state_specs in panel_specs:
        figure = new_figure(498, 140)
        axes = figure.subplots(1, 3)
        for axis, (cell_type, state, state_title) in zip(axes, state_specs):
            subtypes, groups, kruskal_p, radial_p = subtype_groups(states, aoi_type, cell_type, state)
            print(file_name, cell_type, state, "KW", f"{kruskal_p:.2g}", "radial vs vertical", f"{radial_p:.2g}")
            draw_boxes(axis, groups, subtypes, subtype_colours(subtypes))
            axis.set_ylabel(f"{state_title.split(' (')[0]} abundance", fontsize=TITLE_PT)
            title_with_stats(axis, state_title, [f"Kruskal–Wallis $\\it{{P}}$ = {format_p(kruskal_p)}",
                                                  f"radial vs vertical $\\it{{P}}$ = {format_p(radial_p)}"])
        figure.suptitle(row_title, fontsize=TITLE_PT, x=0.01, ha="left")
        save_panel(figure, file_name)


# ================================================================== Fig. 8  cell-type-intrinsic genes
COMPARTMENT_OF = {"Melanoma": "S100", "CAF": "neg_neg", "Endothelial": "neg_neg", "Mural_cell": "neg_neg",
                  "T-CD8": "CD45", "T-CD4": "CD45", "NK": "CD45", "B_cell": "CD45", "Macrophage": "CD45"}
COMPARTMENT_NAME = {"S100": "S100+", "neg_neg": "S100−/CD45−", "CD45": "CD45+"}
# Extra epithelial genes removed when picking heat-map genes (same list as the June selection).
EPITHELIAL_PATTERN = re.compile(r"^(S100A14|JUP|EXPH5|SCNN1A|EPCAM|CDH1|CLDN\d|KRT\d|DSP|PERP|LAD1|TACSTD2|GRHL|PKP\d|CALML\d|SFN|PI3|DMKN|SPRR|LCE|DSG|DSC|FLG|IVL|LOR)$")


def imputed_expression(cell_type):
    # CIBERSORTx high-resolution expression of one cell type in the melanoma AOIs of its compartment.
    compartment = COMPARTMENT_OF[cell_type]
    table = pd.read_csv(f"{HIRES}/CIBERSORTxHiRes_NA_{cell_type}_Window36.txt", sep="\t", index_col=0)
    table.columns = table.columns.astype(str)
    kept_aois = []
    for aoi in table.columns:
        if aoi not in geomx_metadata.index:
            continue
        row = geomx_metadata.loc[aoi]
        if row["Pathophysiology"] == "Melanoma" and row["AOI"] == compartment and row["Path_sub_type"] in SUBTYPES:
            kept_aois.append(aoi)
    expression = table[kept_aois].astype(float)
    expression = expression[(expression == 1).mean(axis=1) < 0.4]     # drop genes not estimated (value 1) in >= 40% of AOIs
    expression = np.log2(expression + 1)
    expression = expression[expression.std(axis=1) > 0]
    subtype_of_aoi = geomx_metadata.loc[kept_aois, "Path_sub_type"]
    aois_by_subtype = {}
    for subtype in SUBTYPES:
        aois = [aoi for aoi in kept_aois if subtype_of_aoi[aoi] == subtype]
        if len(aois) >= 3:
            aois_by_subtype[subtype] = aois
    return compartment, aois_by_subtype, expression


def fig8():
    # Step 1. Fig. 8a: genes with q < 0.05 per cell type (epidermal genes excluded), 150 x 190 pt.
    gene_counts = {}
    for cell_type, compartment in COMPARTMENT_OF.items():
        scan = pd.read_csv(f"{SCAN}/scan_{cell_type}_{compartment}.csv")
        gene_counts[cell_type] = int(((scan["q"] < 0.05) & (~scan["epidermal"].astype(bool))).sum())
    print("8a counts", gene_counts)
    count_order = ["Melanoma", "CAF", "Mural_cell", "Endothelial", "T-CD4", "T-CD8", "NK", "B_cell", "Macrophage"]
    figure = new_figure(150, 190)
    axis = figure.add_subplot()
    y_positions = np.arange(len(count_order))[::-1]
    axis.barh(y_positions, [gene_counts[c] for c in count_order], color=[CELL_COLOUR[c] for c in count_order], height=0.65)
    for y_position, cell_type in zip(y_positions, count_order):
        axis.text(gene_counts[cell_type] + 3, y_position, str(gene_counts[cell_type]), va="center", fontsize=SMALL_PT)
    axis.set_yticks(y_positions)
    axis.set_yticklabels([f"{DISPLAY_NAME[c]} ({COMPARTMENT_NAME[COMPARTMENT_OF[c]]})" for c in count_order], fontsize=SMALL_PT)
    axis.set_xlabel("Genes with q < 0.05", fontsize=TITLE_PT)
    figure.suptitle("Subtype-differential genes", fontsize=TITLE_PT)
    axis.set_xlim(0, max(gene_counts.values()) * 1.25)
    save_panel(figure, "fig8a_counts")

    # Step 2. Fig. 8b: top genes per cell type, z-score of the subtype median (324 x 190 pt).
    genes_shown = [("Melanoma", 6), ("CAF", 4), ("Mural_cell", 4), ("T-CD4", 4), ("Macrophage", 3), ("NK", 3)]
    row_genes = []
    row_cell_types = []
    z_rows = []
    expression_cache = {}
    for cell_type, n_genes in genes_shown:
        compartment, aois_by_subtype, expression = imputed_expression(cell_type)
        expression_cache[cell_type] = (compartment, aois_by_subtype, expression)
        scan = pd.read_csv(f"{SCAN}/scan_{cell_type}_{compartment}.csv")
        scan = scan[(~scan["epidermal"].astype(bool)) & (~scan["gene"].map(lambda g: bool(EPITHELIAL_PATTERN.match(g))))]
        top_genes = scan.sort_values("p").head(n_genes)["gene"].tolist()
        for gene in top_genes:
            medians = np.array([np.median(expression.loc[gene, aois_by_subtype[s]]) for s in SUBTYPES])
            z_rows.append((medians - medians.mean()) / (medians.std() + 1e-9))
            row_genes.append(gene)
            row_cell_types.append(cell_type)
    z_matrix = np.array(z_rows)

    figure = new_figure(324, 190)
    grid = GridSpec(1, 3, figure=figure, width_ratios=[0.6, 12, 0.5], wspace=0.02)
    strip_axis = figure.add_subplot(grid[0, 0])
    heatmap_axis = figure.add_subplot(grid[0, 1])
    colorbar_axis = figure.add_subplot(grid[0, 2])
    image = vector_heatmap(heatmap_axis, z_matrix, "RdBu_r", vmin=-1.5, vmax=1.5)
    heatmap_axis.set_xticks(range(len(SUBTYPES)))
    heatmap_axis.set_xticklabels(SUBTYPES, fontsize=SMALL_PT)
    heatmap_axis.set_yticks(range(len(row_genes)))
    heatmap_axis.set_yticklabels([italic(g) for g in row_genes], fontsize=SMALL_PT)
    heatmap_axis.tick_params(axis="y", length=0, pad=1)
    heatmap_axis.set_title("Top cell-type-intrinsic genes", fontsize=TITLE_PT, pad=14)
    for column, subtype in enumerate(SUBTYPES):
        bar_colour = VERTICAL_COLOUR if VERTICAL_GROWTH[subtype] else RADIAL_COLOUR
        heatmap_axis.add_patch(Rectangle((column - 0.5, -1.45), 1, 0.7, color=bar_colour, clip_on=False))
    heatmap_axis.text(1, -1.1, "Radial growth", color="white", fontsize=SMALL_PT, ha="center", va="center")
    heatmap_axis.text(4, -1.1, "Vertical growth", color="white", fontsize=SMALL_PT, ha="center", va="center")
    block_end = 0
    for _, n_genes in genes_shown:
        block_end = block_end + n_genes
        heatmap_axis.axhline(block_end - 0.5, color="white", lw=1.2)
    strip_axis.set_xlim(0, 1)
    strip_axis.set_ylim(len(row_genes) - 0.5, -0.5)
    strip_axis.axis("off")
    for row_index, cell_type in enumerate(row_cell_types):
        strip_axis.add_patch(Rectangle((0.15, row_index - 0.5), 0.7, 1, color=CELL_COLOUR[cell_type], lw=0))
    strip_axis.legend(handles=[Patch(color=CELL_COLOUR[c], label=DISPLAY_NAME[c]) for c, _ in genes_shown],
                      loc="upper right", bbox_to_anchor=(-0.4, 1.0), frameon=False, fontsize=SMALL_PT,
                      handlelength=0.9, title="Cell type", title_fontsize=TITLE_PT)
    colorbar = keep_colorbar_vector(figure.colorbar(image, cax=colorbar_axis))
    colorbar.set_label("Subtype median (z)", fontsize=TITLE_PT)
    colorbar.ax.tick_params(labelsize=SMALL_PT)
    save_panel(figure, "fig8b_heatmap")

    # Step 3. Fig. 8c: selected candidate genes across subtypes (498 x 250 pt).
    candidates = [("Melanoma", "KIT"), ("Melanoma", "SERPINE2"), ("Melanoma", "AEBP1"), ("Melanoma", "FAT2"),
                  ("Mural_cell", "CFD"), ("T-CD4", "GNG10"), ("Macrophage", "IFI30"), ("NK", "FCGR3A")]
    figure = new_figure(498, 250)
    axes = figure.subplots(2, 4)
    for axis, (cell_type, gene) in zip(axes.flat, candidates):
        if cell_type not in expression_cache:
            expression_cache[cell_type] = imputed_expression(cell_type)
        compartment, aois_by_subtype, expression = expression_cache[cell_type]
        scan = pd.read_csv(f"{SCAN}/scan_{cell_type}_{compartment}.csv").set_index("gene")
        groups = [expression.loc[gene, aois_by_subtype[s]].values for s in SUBTYPES]
        draw_boxes(axis, groups, SUBTYPES, subtype_colours(SUBTYPES), show_n=False)
        print(f"8c {gene} {cell_type}: q={scan.loc[gene, 'q']:.2g} peak {scan.loc[gene, 'peak_subtype']}")
        title_with_stats(axis, f"{italic(gene)} · {DISPLAY_NAME[cell_type]} ({COMPARTMENT_NAME[compartment]})",
                         [f"$\\it{{q}}$ = {format_p(scan.loc[gene, 'q'])}, peak {scan.loc[gene, 'peak_subtype']}"])
        if axis in axes[:, 0]:
            axis.set_ylabel("log2 imputed expression", fontsize=TITLE_PT)
    save_panel(figure, "fig8c_boxplots")


if __name__ == "__main__":
    all_panels = {"6a": fig6a, "6b": fig6b, "6c": fig6c, "6d": fig6d, "6e": fig6e,
                  "7a": fig7a, "7b": fig7b, "7c": fig7c, "7def": fig7def, "8": fig8}
    requested = sys.argv[2:] or list(all_panels.keys())
    for panel_name in requested:
        all_panels[panel_name]()
