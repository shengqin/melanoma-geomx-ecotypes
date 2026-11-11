"""Regenerate matplotlib panels for manuscript Figures 5-8 (vector PDF).
Sources: analysis/ and data/ folders of the GeoMx project. Numbers reproduce the
2026-06-02 follow-up analyses (same AOI sets, same tests)."""
import os, re, sys, numpy as np, pandas as pd, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch, Rectangle
from matplotlib.gridspec import GridSpec
from scipy import stats

BASE = "/Users/Brian/Library/CloudStorage/Box-Box/NanoString GeoMx data analysis"
A = f"{BASE}/analysis"
DISC = f"{A}/2026-01-28_TCGA_SKCM_Discovery_997_4_4"
REC = f"{A}/2026-01-28_GeoMx_recovery_noNorm"
HD = f"{BASE}/data/2026-01-28_melanoma_GeoMx_hires"
OUT = sys.argv[1] if len(sys.argv) > 1 else "panels"
os.makedirs(OUT, exist_ok=True)

plt.rcParams.update({
    "font.family": "Arial", "font.size": 7, "axes.titlesize": 7.5, "axes.labelsize": 7,
    "xtick.labelsize": 6.5, "ytick.labelsize": 6.5, "legend.fontsize": 6.5,
    "axes.spines.top": False, "axes.spines.right": False, "axes.linewidth": 0.6,
    "xtick.major.width": 0.6, "ytick.major.width": 0.6, "xtick.major.size": 2.5, "ytick.major.size": 2.5,
    "pdf.fonttype": 42, "ps.fonttype": 42, "figure.dpi": 100, "savefig.dpi": 300,
})
CELLORD = ["Melanoma", "T-CD8", "T-CD4", "NK", "B_cell", "Macrophage", "Endothelial", "CAF", "Mural_cell"]
NICE = {"Melanoma": "Melanoma", "T-CD8": "T-CD8", "T-CD4": "T-CD4", "NK": "NK", "B_cell": "B cell",
        "Macrophage": "Macrophage", "Endothelial": "Endothelial", "CAF": "CAF", "Mural_cell": "Mural cell"}
PAL = {"Melanoma": "#7B4A23", "T-CD8": "#2C6FB0", "T-CD4": "#7FB2E0", "NK": "#9C27B0", "B_cell": "#3FA34D",
       "Macrophage": "#E08A2B", "Endothelial": "#D7263D", "CAF": "#C85A54", "Mural_cell": "#6D597A"}
SUBS = ["LMM", "SITU", "SSM", "NOD", "AMEL", "MET"]
RAD = {"LMM": 0, "SITU": 0, "SSM": 0, "NOD": 1, "AMEL": 1, "MET": 1}
C_RAD, C_VERT = "#3E7C7C", "#C85A54"
E1C, E2C = "#D7263D", "#F2C14E"




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


def save(fig, name):
    p = f"{OUT}/{name}.pdf"
    fig.savefig(p)
    fig.savefig(p.replace(".pdf", ".png"), dpi=150)
    plt.close(fig)
    print("saved", name)


def fmt_p(p):
    if p < 1e-3:
        m, e = f"{p:.0e}".split("e")
        return f"{m}×10$^{{{int(e)}}}$"
    return f"{p:.3f}" if p < 0.01 else f"{p:.2f}"


def boxes(ax, data, labels, colors, ylabel, title=None, show_n=True, point_size=6):
    bp = ax.boxplot(data, patch_artist=True, widths=0.62, showfliers=False,
                    medianprops=dict(color="black", lw=0.9), whiskerprops=dict(lw=0.6), capprops=dict(lw=0.6),
                    boxprops=dict(lw=0.6))
    for b, c in zip(bp["boxes"], colors):
        b.set_facecolor(c); b.set_alpha(0.75)
    rng = np.random.default_rng(1)
    for i, d in enumerate(data):
        ax.scatter(rng.normal(i + 1, 0.07, len(d)), d, s=point_size, c="k", alpha=0.55, lw=0, zorder=3)
    xl = [f"{l}\n(n={len(d)})" if show_n else l for l, d in zip(labels, data)]
    ax.set_xticks(range(1, len(labels) + 1)); ax.set_xticklabels(xl)
    ax.set_ylabel(ylabel)
    if title: ax.set_title(title, fontsize=7.5)
    return bp


# ------------------------------------------------------------------ metadata
meta = pd.read_csv(f"{BASE}/data/2025/2025-01-27_GeoMxmetadata.txt", sep="\t")
meta["ID"] = meta["ID"].astype(str); meta = meta.set_index("ID")
meta["AOI"] = meta["Segment_Tags"].replace({"neg/neg": "neg_neg"})
asg = pd.read_csv(f"{REC}/Ecotypes/ecotype_assignment.txt", sep="\t", index_col=0); asg.index = asg.index.astype(str)
asg["AOI"] = asg["Segment_Tags"].replace({"neg/neg": "neg_neg"})


def load_states(folder, cell):
    d = pd.read_csv(f"{folder}/{cell}/state_abundances.txt", sep="\t", index_col=0)
    d.columns = d.columns.astype(str)
    return d.T


# ================================================================== FIG 5b full signature-matrix chessboard
# Typography standard for Fig. 5 (author request, 30 Oct 2026): Arial, 7 pt titles, 6 pt ticks and small text.
# Fig. 5 panels are drawn at their exact printed size (in points) so compose.py places them unscaled.
TITLE_PT = 7
SMALL_PT = 6
POINTS_PER_INCH = 72


def fig5b():
    # Step 1. Read the signature matrix and use the same cell-type names and order as the rest of Fig. 5.
    signature = pd.read_csv(f"{BASE}/data/2026-01-20_melanoma_sigmtx_400_800.txt", sep="\t", index_col=0)
    column_renames = {"B-cell": "B_cell", "Mural\u2212cell": "Mural_cell", "Mural-cell": "Mural_cell"}
    signature.columns = [column_renames.get(name, name) for name in signature.columns]
    signature = signature[CELLORD].astype(float)

    # Step 2. Row z-score of log2 expression, so each gene shows where it peaks.
    log_expression = np.log2(signature + 1)
    row_means = log_expression.mean(axis=1)
    row_sds = log_expression.std(axis=1) + 1e-9          # small constant avoids division by zero
    z_scores = log_expression.sub(row_means, axis=0).div(row_sds, axis=0)

    # Step 3. Group genes by the cell type where they peak, strongest first within each group.
    peak_cell_type = z_scores.idxmax(axis=1)
    gene_order = []
    block_ends = []
    for cell_type in CELLORD:
        genes_in_block = z_scores.index[peak_cell_type == cell_type]
        genes_in_block = z_scores.loc[genes_in_block, cell_type].sort_values(ascending=False).index.tolist()
        gene_order = gene_order + genes_in_block
        block_ends.append(len(gene_order))
    heatmap_values = z_scores.loc[gene_order].values
    print("5b genes:", len(gene_order))

    # Step 4. Draw at the printed size: 140 x 200 pt (same height as the UMAPs in row 1).
    figure = plt.figure(figsize=(140 / POINTS_PER_INCH, 200 / POINTS_PER_INCH), layout="constrained")
    grid = GridSpec(1, 2, figure=figure, width_ratios=[16, 1], wspace=0.04)
    heatmap_axis = figure.add_subplot(grid[0, 0])
    colorbar_axis = figure.add_subplot(grid[0, 1])
    image = vector_heatmap(heatmap_axis, heatmap_values, "RdBu_r", vmin=-2, vmax=2)
    for block_end in block_ends[:-1]:
        heatmap_axis.axhline(block_end - 0.5, color="k", lw=0.4)

    # Step 5. Labels: 7 pt titles, 6 pt ticks.
    heatmap_axis.set_title("Signature matrix", fontsize=TITLE_PT)
    heatmap_axis.set_ylabel(f"Signature genes (n = {len(gene_order):,})", fontsize=TITLE_PT)
    heatmap_axis.set_yticks([])
    heatmap_axis.set_xticks(range(len(CELLORD)))
    heatmap_axis.set_xticklabels([NICE[c] for c in CELLORD], rotation=90, ha="center", fontsize=SMALL_PT)   # vertical: 45 degrees overlapped
    colorbar = keep_colorbar_vector(figure.colorbar(image, cax=colorbar_axis))
    colorbar.set_label("Row z-score (log2 expression)", fontsize=TITLE_PT)
    colorbar.ax.tick_params(labelsize=SMALL_PT)
    save(figure, "fig5b_sigmatrix_full")


# ================================================================== FIG 5c pseudobulk scatter
def fig5c():
    d = f"{BASE}/data/2026-01-28_melanoma_pseudobulk"
    cells = ["Melanoma", "T-CD8", "T-CD4", "NK", "B_cell", "Macrophage", "Endothelial", "CAF", "Mural_cell"]
    fig, axes = plt.subplots(3, 3, figsize=(4.05, 3.5), layout="constrained")
    for ax, cell in zip(axes.flat, cells):
        r = pd.read_csv(f"{d}/2026-01-28_melanoma_sigmtx_400_800_{cell}_pseudo/CIBERSORTx_Results.txt", sep="\t")
        r = r[r["Mixture"] != "sample_control"].copy()
        r["truth"] = r["Mixture"].str[-3:].astype(float) / 1000
        est = r[cell].astype(float).values; tr = r["truth"].values
        rho, p = stats.spearmanr(tr, est)
        # Lin's concordance
        mx, my = tr.mean(), est.mean(); vx, vy = tr.var(), est.var(); cov = np.mean((tr - mx) * (est - my))
        ccc = 2 * cov / (vx + vy + (mx - my) ** 2)
        lim = max(tr.max(), est.max()) * 1.08
        ax.plot([0, lim], [0, lim], color="grey", lw=0.6, ls="--")
        ax.scatter(tr, est, s=9, color=PAL[cell], edgecolor="k", lw=0.2, zorder=3)
        ax.set_title(NICE[cell], fontsize=7, pad=2)
        ax.text(0.04, 0.96, f"ρ={rho:.2f}\nCCC={ccc:.2f}", transform=ax.transAxes, va="top", fontsize=5.8)
        ax.set_xlim(0, lim); ax.set_ylim(0, lim)
        ax.tick_params(labelsize=5.5, pad=1.5)
        ax.locator_params(nbins=4)
    for ax in axes[2, :]: ax.set_xlabel("Spiked-in fraction (truth)", fontsize=6.3)
    for ax in axes[:, 0]: ax.set_ylabel("Estimated fraction", fontsize=6.3)
    save(fig, "fig5c_pseudobulk_scatter")


# ================================================================== FIG 5c pseudobulk Spearman rho bar plot
def fig5c_rho():
    # Step 1. Spearman correlation between estimated and spiked-in fraction, one value per cell type.
    pseudobulk_dir = f"{BASE}/data/2026-01-28_melanoma_pseudobulk"
    cell_types = ["Melanoma", "T-CD8", "T-CD4", "NK", "B_cell", "Macrophage", "Endothelial", "CAF", "Mural_cell"]
    spearman_values = []
    for cell_type in cell_types:
        results = pd.read_csv(f"{pseudobulk_dir}/2026-01-28_melanoma_sigmtx_400_800_{cell_type}_pseudo/CIBERSORTx_Results.txt", sep="\t")
        results = results[results["Mixture"] != "sample_control"].copy()     # the control has no spike-in
        results["truth"] = results["Mixture"].str[-3:].astype(float) / 1000    # mixture name ends in the spiked fraction x 1000
        rho, _ = stats.spearmanr(results["truth"].values, results[cell_type].astype(float).values)
        spearman_values.append(rho)
    print("5c rho:", dict(zip(cell_types, [round(v, 3) for v in spearman_values])))

    # Step 2. Bar plot at the printed size: 245 x 150 pt.
    figure, axis = plt.subplots(figsize=(245 / POINTS_PER_INCH, 150 / POINTS_PER_INCH), layout="constrained")
    positions = np.arange(len(cell_types))
    axis.bar(positions, spearman_values, color=[PAL[c] for c in cell_types], width=0.65, edgecolor="k", lw=0.4)
    for position, value in zip(positions, spearman_values):
        axis.text(position, value + 0.015, f"{value:.2f}", ha="center", va="bottom", fontsize=SMALL_PT)

    # Step 3. Labels: 7 pt titles, 6 pt ticks.
    axis.set_title("Pseudobulk benchmark (21 mixtures per cell type)", fontsize=TITLE_PT)
    axis.set_ylabel("Spearman \u03c1 (estimated vs\nspiked-in fraction)", fontsize=TITLE_PT)
    axis.set_xlabel("Cell type", fontsize=TITLE_PT)
    axis.set_xticks(positions)
    axis.set_xticklabels([NICE[c] for c in cell_types], rotation=40, ha="right", fontsize=SMALL_PT)
    axis.set_ylim(0, 1.1)
    axis.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    axis.tick_params(axis="y", labelsize=SMALL_PT)
    save(figure, "fig5c_pseudobulk_rho")


# ================================================================== FIG 6a composition
def fig6a():
    fr = pd.read_csv(f"{BASE}/data/2026-01-28_TCGA_SKCM_fraction.txt", sep="\t", index_col=0)
    F = fr[CELLORD].astype(float); F = F.div(F.sum(1), axis=0).sort_values("Melanoma")
    fig, ax = plt.subplots(figsize=(7.08, 1.95), layout="constrained")
    x = np.arange(F.shape[0]); bottom = np.zeros(F.shape[0])
    for ct in CELLORD:
        ax.bar(x, F[ct].values, bottom=bottom, width=1.0, color=PAL[ct], label=NICE[ct], linewidth=0)
        bottom += F[ct].values
    ax.set_xlim(-0.5, F.shape[0] - 0.5); ax.set_ylim(0, 1)
    ax.set_xlabel(f"TCGA-SKCM tumours ranked by melanoma fraction (n={F.shape[0]})")
    ax.set_ylabel("Cell-type fraction\n(CIBERSORTx)")
    ax.legend(ncol=9, loc="upper center", bbox_to_anchor=(0.5, -0.32), frameon=False, handlelength=1, columnspacing=1.2)
    save(fig, "fig6a_tcga_composition")


# ================================================================== FIG 6b melanoma state heatmap
def fig6b():
    H = pd.read_csv(f"{DISC}/Melanoma/heatmap_data.txt", sep="\t", index_col=0); H.columns = H.columns.astype(str)
    gi = pd.read_csv(f"{DISC}/Melanoma/gene_info.txt", sep="\t", index_col=0)
    ann = pd.read_csv(f"{DISC}/Melanoma/heatmap_top_ann.txt", sep="\t", index_col=0); ann.index = ann.index.astype(str)
    st = load_states(DISC, "Melanoma")
    samples = [s for s in H.columns if s in st.index]
    order = sorted(samples, key=lambda s: (-(ann.loc[s, "State"] == "S01"), -st.loc[s, "S01"]))
    g1 = [g for g in H.index if gi.loc[g, "State"] == "S01"]; g2 = [g for g in H.index if gi.loc[g, "State"] == "S02"]
    genes = g1 + g2
    M = H.loc[genes, order].values
    M = np.clip(M, -2.5, 2.5)
    fig = plt.figure(figsize=(3.5, 3.2), layout="constrained")
    gs = GridSpec(4, 2, figure=fig, height_ratios=[0.5, 0.5, 0.5, 12], width_ratios=[30, 3.6], hspace=0.04, wspace=0.0)
    n = len(order)
    def bar(row, vals, cmap, label):
        ax = fig.add_subplot(gs[row, 0]); ax.imshow(vals[None, :], aspect="auto", cmap=cmap, interpolation="nearest")
        ax.set_xticks([]); ax.set_yticks([]); [s.set_visible(False) for s in ax.spines.values()]
        ax.text(-0.005, 0.5, label, transform=ax.transAxes, ha="right", va="center", fontsize=6)
        return ax
    from matplotlib.colors import ListedColormap
    state_v = np.array([0 if ann.loc[s, "State"] == "S01" else 1 for s in order])
    bar(0, state_v, ListedColormap(["#E8336D", "#4BA3A3"]), "State")
    typ = ann.loc[order, "SAMPLE_TYPE"].astype(str)
    typ_v = np.array([0 if t == "Primary" else (1 if t == "Metastasis" else 2) for t in typ])
    bar(1, typ_v, ListedColormap(["#F2C14E", "#8C8C8C", "#EEEEEE"]), "Sample type")
    s01 = st.loc[order, "S01"].values
    bar(2, s01, "viridis", "S01 abund.")
    ax = fig.add_subplot(gs[3, 0])
    im = ax.imshow(M, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5, interpolation="nearest")
    ax.axhline(len(g1) - 0.5, color="k", lw=0.8)
    ax.axvline(int((state_v == 0).sum()) - 0.5, color="k", lw=0.8)
    ax.set_xticks([]); ax.set_yticks([])
    ax.set_xlabel(f"TCGA-SKCM tumours (n={n}), ordered by state and S01 abundance", fontsize=6.3)
    ax.text(-0.02, 1 - len(g1) / 2 / len(genes), f"S01\npigmented\n({len(g1)} genes)", transform=ax.transAxes, ha="right", va="center", fontsize=6, color="#E8336D", fontweight="bold")
    ax.text(-0.02, (len(g2) / 2) / len(genes), f"S02\nde-diff.\n({len(g2)} genes)", transform=ax.transAxes, ha="right", va="center", fontsize=6, color="#2E7D7D", fontweight="bold")
    # marker labels on right: leader lines drawn from the heatmap edge into the label column
    want1 = ["OCA2", "PMEL", "GPR143", "TRPM1", "TSPAN10", "MFSD12", "MLANA", "TYR", "MITF"]
    want2 = ["AXL", "NGFR", "SOX9", "WNT5A", "TGFBI", "SERPINE1", "IL6ST", "AHR", "LEF1", "CD68", "AEBP1", "SERPINE2"]
    axl = fig.add_subplot(gs[3, 1]); axl.set_xlim(0, 1); axl.set_ylim(len(genes) - 0.5, -0.5); axl.axis("off")
    lab = sorted((genes.index(g), g) for g in want1 + want2 if g in genes)
    ys = np.array([i for i, _ in lab], float); min_gap = len(genes) / 26
    for k in range(1, len(ys)):
        if ys[k] - ys[k - 1] < min_gap: ys[k] = ys[k - 1] + min_gap
    for (i, g), y in zip(lab, ys):
        ax.annotate(g, xy=(n - 0.5, i), xycoords="data", xytext=(1.06, 1 - (y + 0.5) / len(genes)), textcoords="axes fraction",
                    fontsize=5.5, style="italic", va="center", ha="left", annotation_clip=False,
                    arrowprops=dict(arrowstyle="-", lw=0.4, color="grey", shrinkA=0, shrinkB=0))
    cb = fig.colorbar(im, ax=ax, location="bottom", fraction=0.03, pad=0.04, shrink=0.28, anchor=(1.0, 0))
    cb.set_label("expression (z)", fontsize=6); cb.ax.tick_params(labelsize=5.5)
    fig.legend(handles=[Patch(color="#E8336D", label="S01"), Patch(color="#4BA3A3", label="S02"),
                        Patch(color="#F2C14E", label="Primary"), Patch(color="#8C8C8C", label="Metastasis"), Patch(color="#EEEEEE", label="not annotated")],
               loc="outside lower left", ncol=5, frameon=False, fontsize=5.8, handlelength=1, columnspacing=0.9, handletextpad=0.4)
    save(fig, "fig6b_tcga_melanoma_states")
    return order


# ================================================================== FIG 6c S01 vs signatures
PIG = ["PMEL", "MLANA", "TYR", "TYRP1", "DCT", "MITF", "OCA2", "GPR143", "TRPM1", "SLC45A2", "MFSD12", "TSPAN10"]
IFN = ["STAT1", "IRF1", "CXCL9", "CXCL10", "CXCL11", "GBP1", "IDO1", "HLA-DRA", "CD8A", "GZMB", "PRF1", "CD274", "TAP1"]


def fig6c():
    want = set(PIG + IFN); rows = []; hdr = None
    with open(f"{BASE}/data/2025/2025-01-29_SKCM_coding_symbol_TPM.txt") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            g = line.split("\t", 1)[0]
            if g in want: rows.append(line.rstrip("\n").split("\t"))
    tpm = pd.DataFrame([r[1:] for r in rows], index=[r[0] for r in rows], columns=hdr[1:]).astype(float)
    st = load_states(DISC, "Melanoma")
    common = [s for s in st.index if s in tpm.columns]
    tpm = np.log2(tpm[common] + 1)
    zz = tpm.sub(tpm.mean(1), axis=0).div(tpm.std(1) + 1e-9, axis=0)
    pig = zz.loc[[g for g in PIG if g in zz.index]].mean(0); ifn = zz.loc[[g for g in IFN if g in zz.index]].mean(0)
    s01 = st.loc[common, "S01"]
    fig, axes = plt.subplots(2, 1, figsize=(2.55, 3.0), layout="constrained")
    for ax, sig, lab in [(axes[0], pig, "Pigmentation/MITF\nsignature (z)"), (axes[1], ifn, "IFN-γ/inflammation\nsignature (z)")]:
        r, p = stats.spearmanr(s01, sig)
        ax.scatter(s01, sig, c=s01, cmap="RdBu", vmin=0.15, vmax=0.85, s=7, edgecolor="k", lw=0.15)
        m, b = np.polyfit(s01, sig, 1); xs = np.linspace(0, 1, 5); ax.plot(xs, m * xs + b, color="grey", ls="--", lw=0.7)
        ax.set_ylabel(lab); ax.set_xlim(-0.02, 1.02)
        if r > 0:
            ax.text(0.03, 0.97, f"Spearman ρ = {r:.2f}\nP = {fmt_p(p)}", transform=ax.transAxes, fontsize=6.2, va="top", ha="left")
        else:
            ax.text(0.97, 0.97, f"Spearman ρ = {r:.2f}\nP = {fmt_p(p)}", transform=ax.transAxes, fontsize=6.2, va="top", ha="right")
    axes[1].set_xlabel(f"Melanoma-S01 (pigmented) abundance\nTCGA-SKCM (n={len(common)})")
    axes[0].tick_params(labelbottom=False)
    save(fig, "fig6c_s01_vs_signatures")


# ================================================================== FIG 6e ecotype associations
def fig6e():
    d = pd.read_csv(f"{DISC}/Ecotypes/ecotype_assignment.txt", sep="\t", index_col=0)
    d["ANEUPLOIDY_SCORE"] = pd.to_numeric(d["ANEUPLOIDY_SCORE"], errors="coerce")
    d["TMB_NONSYNONYMOUS"] = pd.to_numeric(d["TMB_NONSYNONYMOUS"], errors="coerce")
    fig, axes = plt.subplots(1, 3, figsize=(3.0, 1.8), layout="constrained", width_ratios=[1.3, 1, 1])
    # sample type
    ax = axes[0]
    dd = d[d.SAMPLE_TYPE.isin(["Primary", "Metastasis"])]
    ct = pd.crosstab(dd.Ecotype, dd.SAMPLE_TYPE)[["Primary", "Metastasis"]]
    OR, p = stats.fisher_exact(ct.values)  # rows E1,E2 ; cols Primary, Metastasis
    pct = ct.div(ct.sum(1), axis=0) * 100
    ax.bar([0, 1], pct["Metastasis"], color="#8C8C8C", width=0.65, label="Metastasis")
    ax.bar([0, 1], pct["Primary"], bottom=pct["Metastasis"], color="#F2C14E", width=0.65, label="Primary")
    ax.set_xticks([0, 1]); ax.set_xticklabels([f"E1\nn={ct.loc['E1'].sum()}", f"E2\nn={ct.loc['E2'].sum()}"], fontsize=6)
    ax.set_ylabel("% of tumours"); ax.set_ylim(0, 100)
    ax.set_title(f"Sample type\nFisher P = {fmt_p(p)}", fontsize=6.5)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.32), frameon=False, ncol=1, fontsize=5.5, handlelength=0.9)
    print(f"6e Fisher: n={len(dd)} OR(E1 vs E2, met odds)={1/OR:.2f} p={p:.2e}", ct.to_dict())
    for ax, var, lab, log in [(axes[1], "ANEUPLOIDY_SCORE", "Aneuploidy score", False), (axes[2], "TMB_NONSYNONYMOUS", "TMB (mut/Mb)", True)]:
        a = d[d.Ecotype == "E1"][var].dropna(); b = d[d.Ecotype == "E2"][var].dropna()
        p = stats.mannwhitneyu(a, b).pvalue
        boxes(ax, [a.values, b.values], ["E1", "E2"], [E1C, E2C], lab, show_n=False, point_size=3)
        ax.set_title(f"{lab.split(' (')[0]}\nP = {fmt_p(p)} (n.s.)", fontsize=6.5)
        if log: ax.set_yscale("log")
        ax.set_ylabel(lab)
    save(fig, "fig6e_tcga_associations")


# ================================================================== FIG 7b-d GeoMx states by subtype
def fig7bcd():
    S = {c: load_states(REC, c) for c in ["Melanoma", "Endothelial", "CAF", "Mural_cell", "B_cell", "T.CD8", "NK", "T.CD4", "Macrophage"]}
    mel = asg[asg.Pathophysiology == "Melanoma"]

    def sub_data(aoi, cell, state):
        idx = mel[mel.AOI == aoi].index
        ser = S[cell][state].reindex(idx)
        df = pd.DataFrame({"v": ser, "sub": mel.loc[idx, "Path_sub_type"]}).dropna()
        order = [s for s in SUBS if (df["sub"] == s).sum() >= 3]
        data = [df[df["sub"] == s]["v"].values for s in order]
        H, p = stats.kruskal(*data)
        df["grp"] = df["sub"].map(RAD)
        pmw = stats.mannwhitneyu(df[df.grp == 0]["v"], df[df.grp == 1]["v"]).pvalue
        return order, data, p, pmw

    # 7b
    order, data, p, pmw = fig7b_data = sub_data("S100", "Melanoma", "S01")
    fig, ax = plt.subplots(figsize=(3.3, 2.6), layout="constrained")
    boxes(ax, data, order, [C_RAD if RAD[s] == 0 else C_VERT for s in order], "Melanoma-S01 (pigmented)\nstate abundance")
    ax.axhline(0.5, color="grey", ls="--", lw=0.6)
    ax.set_title(f"S100+ tumour AOIs   Kruskal–Wallis P = {fmt_p(p)}", fontsize=7)
    ax.set_xlabel("Melanoma histologic subtype")
    ax.legend(handles=[Patch(color=C_RAD, label="radial growth"), Patch(color=C_VERT, label="vertical growth / metastasis")],
              loc="lower left", frameon=False, fontsize=6, handlelength=1)
    print("7b medians:", {s: round(float(np.median(d)), 2) for s, d in zip(order, data)}, "P", p)
    save(fig, "fig7b_S01_by_subtype")

    def triple(name, aoi, spec, ylab_prefix, title):
        fig, axes = plt.subplots(1, 3, figsize=(7.08, 2.05), layout="constrained")
        for ax, (cell, state, lab) in zip(axes, spec):
            order, data, p, pmw = sub_data(aoi, cell, state)
            boxes(ax, data, order, [C_RAD if RAD[s] == 0 else C_VERT for s in order], f"{lab.split(' (')[0]} abundance")
            ax.set_title(f"{lab}\nKruskal–Wallis P = {fmt_p(p)}\nradial vs vertical P = {fmt_p(pmw)}", fontsize=6.3)
            print(name, cell, state, {s: round(float(np.median(d)), 2) for s, d in zip(order, data)}, "KW", p, "MW", pmw)
        fig.suptitle(title, fontsize=7, fontweight="bold", x=0.01, ha="left")
        save(fig, name)

    triple("fig7c_stromal_states", "neg_neg",
           [("Endothelial", "S01", "Endothelial-S01 (quiescent vessel)"), ("CAF", "S01", "CAF-S01 (myofibroblast/desmoplastic)"),
            ("Mural_cell", "S01", "Mural-S01 (mature pericyte)")], "", "Stromal (S100−/CD45−) AOIs")
    triple("fig7d_immune_states", "CD45",
           [("B_cell", "S02", "B cell-S02"), ("T.CD8", "S01", "CD8 T-S01 (memory)"), ("NK", "S01", "NK-S01")], "",
           "Immune (CD45+) AOIs — trend level, none significant after FDR")


# ================================================================== FIG 7b ecotype composition by segment (melanoma AOIs)
def fig7b_eco():
    mel = asg[asg.Pathophysiology == "Melanoma"]
    segs = [("S100", "Tumour\nS100+"), ("CD45", "Immune\nCD45+"), ("neg_neg", "Stromal\nS100−\nCD45−")]
    ct = pd.crosstab(mel.AOI, mel.Ecotype).reindex([s for s, _ in segs])
    pct = ct.div(ct.sum(1), axis=0) * 100
    print("7b ecotype by segment (melanoma AOIs):", ct.to_dict(), pct.round(1).to_dict())
    fig, ax = plt.subplots(figsize=(1.7, 2.6), layout="constrained")
    x = np.arange(len(segs))
    ax.bar(x, pct["E2"], color=E2C, width=0.65, label="E2", edgecolor="k", lw=0.4)
    ax.bar(x, pct["E1"], bottom=pct["E2"], color=E1C, width=0.65, label="E1", edgecolor="k", lw=0.4)
    for xi, seg in zip(x, [s for s, _ in segs]):
        ax.text(xi, pct.loc[seg, "E2"] / 2, f"{pct.loc[seg, 'E2']:.1f}%", ha="center", va="center", fontsize=6)
        ax.text(xi, 101, f"n={int(ct.loc[seg].sum())}", ha="center", va="bottom", fontsize=5.8)
    ax.set_xticks(x); ax.set_xticklabels([l for _, l in segs], fontsize=6)
    ax.set_ylim(0, 108); ax.set_yticks([0, 25, 50, 75, 100]); ax.set_ylabel("% of melanoma AOIs")
    ax.set_title("Ecotype by segment", fontsize=7)
    ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.30), ncol=2, frameon=False, fontsize=6, handlelength=1)
    save(fig, "fig7b_ecotype_by_segment")


# ================================================================== FIG 7e niche coupling
def fig7e():
    F = f"{A}/2026-06-02_ecotype_recovery_followup"
    roi = pd.read_csv(f"{F}/A2_roi_ecotype_by_AOI.csv", index_col=0)
    pw = pd.read_csv(f"{F}/A2_pairwise_concordance.csv")
    fig, axes = plt.subplots(1, 2, figsize=(3.4, 2.35), layout="constrained", width_ratios=[1.1, 1])
    ax = axes[0]
    cols = ["S100", "CD45", "neg_neg"]
    M = (roi[cols] == "E1").astype(int).values
    o = np.lexsort((M[:, 2], M[:, 1], M[:, 0]))[::-1]
    from matplotlib.colors import ListedColormap
    ax.imshow(M[o], aspect="auto", cmap=ListedColormap([E2C, E1C]), interpolation="nearest")
    ax.set_xticks(range(3)); ax.set_xticklabels(["S100+", "CD45+", "S100−/\nCD45−"])
    ax.set_yticks([]); ax.set_ylabel(f"Regions with all 3 AOIs (n={len(roi)})")
    conc = int(roi["concordant"].sum())
    ax.set_title(f"Ecotype per AOI within region\n{conc}/{len(roi)} concordant ({100*conc/len(roi):.0f}%; expected 25%)", fontsize=6.2)
    ax.legend(handles=[Patch(color=E1C, label="E1"), Patch(color=E2C, label="E2")], loc="upper left", bbox_to_anchor=(1.0, 1.0), frameon=False, fontsize=6, handlelength=0.9)
    for s in ax.spines.values(): s.set_visible(True)
    ax = axes[1]
    ax.bar(range(3), pw["agreement_pct"], color="#8E7CC3", width=0.6)
    ax.axhline(50, color="grey", ls="--", lw=0.6)
    for i, (a, k) in enumerate(zip(pw["agreement_pct"], pw["cohen_kappa"])):
        ax.text(i, a + 2, f"κ={k:.2f}", ha="center", fontsize=5.8)
    ax.set_xticks(range(3)); ax.set_xticklabels(["S100 vs\nCD45", "S100 vs\nstromal", "CD45 vs\nstromal"])
    ax.set_ylim(0, 100); ax.set_ylabel("Pairwise agreement (%)")
    ax.set_title("Pairwise ecotype agreement\n(dashed: chance)", fontsize=6.3)
    save(fig, "fig7e_niche_coupling")


# ================================================================== FIG 8 hires candidates
PAIR = {"Melanoma": "S100", "CAF": "neg_neg", "Endothelial": "neg_neg", "Mural_cell": "neg_neg",
        "T-CD8": "CD45", "T-CD4": "CD45", "NK": "CD45", "B_cell": "CD45", "Macrophage": "CD45"}
EXTRA = re.compile(r"^(S100A14|JUP|EXPH5|SCNN1A|EPCAM|CDH1|CLDN\d|KRT\d|DSP|PERP|LAD1|TACSTD2|GRHL|PKP\d|CALML\d|SFN|PI3|DMKN|SPRR|LCE|DSG|DSC|FLG|IVL|LOR)$")
SCAN = f"{A}/2026-06-02_hires_subtype_genes"


def hires_matrix(cell):
    comp = PAIR[cell]
    df = pd.read_csv(f"{HD}/CIBERSORTxHiRes_NA_{cell}_Window36.txt", sep="\t", index_col=0); df.columns = df.columns.astype(str)
    keep = [c for c in df.columns if c in meta.index and meta.loc[c, "Pathophysiology"] == "Melanoma"
            and meta.loc[c, "AOI"] == comp and meta.loc[c, "Path_sub_type"] in SUBS]
    sub = meta.loc[keep, "Path_sub_type"]
    X = df[keep].astype(float); X = X[(X == 1).mean(1) < 0.4]; X = np.log2(X + 1); X = X[X.std(1) > 0]
    grp = {s: [c for c in keep if sub[c] == s] for s in SUBS if (sub == s).sum() >= 3}
    return comp, grp, X


def fig8():
    counts = {}
    for cell in PAIR:
        R = pd.read_csv(f"{SCAN}/scan_{cell}_{PAIR[cell]}.csv")
        counts[cell] = int(((R.q < 0.05) & (~R.epidermal.astype(bool))).sum())
    print("8a counts", counts)
    # ---- a: counts bar
    fig, ax = plt.subplots(figsize=(2.05, 2.55), layout="constrained")
    order = ["Melanoma", "CAF", "Mural_cell", "Endothelial", "T-CD4", "T-CD8", "NK", "B_cell", "Macrophage"]
    comp_lab = {"S100": "S100+", "neg_neg": "S100−/CD45−", "CD45": "CD45+"}
    y = np.arange(len(order))[::-1]
    ax.barh(y, [counts[c] for c in order], color=[PAL[c] for c in order], height=0.65)
    for yi, c in zip(y, order):
        ax.text(counts[c] + 3, yi, str(counts[c]), va="center", fontsize=6)
    ax.set_yticks(y); ax.set_yticklabels([f"{NICE[c]} ({comp_lab[PAIR[c]]})" for c in order], fontsize=6)
    ax.set_xlabel("Genes with q < 0.05\n(cell-type-intrinsic)", fontsize=6.5)
    ax.set_xlim(0, max(counts.values()) * 1.22)
    save(fig, "fig8a_counts")
    # ---- b: heatmap
    SHOW = [("Melanoma", 6), ("CAF", 4), ("Mural_cell", 4), ("T-CD4", 4), ("Macrophage", 3), ("NK", 3)]
    rows, rowcell, zmat, mats = [], [], [], {}
    for cell, n in SHOW:
        comp, grp, X = hires_matrix(cell); mats[cell] = (comp, grp, X)
        R = pd.read_csv(f"{SCAN}/scan_{cell}_{comp}.csv")
        R = R[(~R.epidermal.astype(bool)) & (~R.gene.map(lambda g: bool(EXTRA.match(g))))]
        top = R.sort_values("p").head(n)["gene"].tolist()
        for g in top:
            m = np.array([np.median(X.loc[g, grp[s]]) for s in SUBS])
            zmat.append((m - m.mean()) / (m.std() + 1e-9)); rows.append(g); rowcell.append(cell)
    Z = np.array(zmat)
    fig = plt.figure(figsize=(4.5, 2.55), layout="constrained")
    gs = GridSpec(1, 3, figure=fig, width_ratios=[0.5, 12, 0.5], wspace=0.02)
    axs = fig.add_subplot(gs[0, 0]); ax = fig.add_subplot(gs[0, 1]); axc = fig.add_subplot(gs[0, 2])
    im = ax.imshow(Z, aspect="auto", cmap="RdBu_r", vmin=-1.5, vmax=1.5, interpolation="nearest")
    ax.set_xticks(range(len(SUBS))); ax.set_xticklabels(SUBS)
    ax.set_yticks(range(len(rows))); ax.set_yticklabels(rows, fontsize=5.8, style="italic")
    ax.tick_params(axis="y", length=0, pad=1)
    for j, s in enumerate(SUBS):
        ax.add_patch(Rectangle((j - 0.5, -1.45), 1, 0.7, color=C_RAD if RAD[s] == 0 else C_VERT, clip_on=False))
    ax.text(1, -1.75, "radial growth", color=C_RAD, fontsize=6, ha="center")
    ax.text(4, -1.75, "vertical growth", color=C_VERT, fontsize=6, ha="center")
    b = 0
    for cell, n in SHOW:
        b += n; ax.axhline(b - 0.5, color="white", lw=1.2)
    axs.set_xlim(0, 1); axs.set_ylim(len(rows) - 0.5, -0.5); axs.axis("off")
    for i, c in enumerate(rowcell): axs.add_patch(Rectangle((0.15, i - 0.5), 0.7, 1, color=PAL[c], lw=0))
    axs.legend(handles=[Patch(color=PAL[c], label=NICE[c]) for c, _ in SHOW], loc="upper right", bbox_to_anchor=(-0.4, 1.0), frameon=False, fontsize=5.8, handlelength=0.9, title="Cell type", title_fontsize=6)
    cb = fig.colorbar(im, cax=axc); cb.set_label("z of subtype median", fontsize=6); cb.ax.tick_params(labelsize=5.5)
    save(fig, "fig8b_heatmap")
    # ---- c: boxplots
    PANELS = [("Melanoma", "KIT"), ("Melanoma", "SERPINE2"), ("Melanoma", "AEBP1"), ("Melanoma", "FAT2"),
              ("Mural_cell", "CFD"), ("T-CD4", "GNG10"), ("Macrophage", "IFI30"), ("NK", "FCGR3A")]
    fig, axes = plt.subplots(2, 4, figsize=(7.08, 3.5), layout="constrained")
    for ax, (cell, g) in zip(axes.flat, PANELS):
        if cell not in mats: mats[cell] = hires_matrix(cell)
        comp, grp, X = mats[cell]
        R = pd.read_csv(f"{SCAN}/scan_{cell}_{comp}.csv").set_index("gene")
        data = [X.loc[g, grp[s]].values for s in SUBS]
        boxes(ax, data, SUBS, [C_RAD if RAD[s] == 0 else C_VERT for s in SUBS], "$\\log_2$ imputed expression", show_n=False, point_size=5)
        ax.set_title(f"$\\it{{{g}}}$ · {NICE[cell]} ({comp_lab[comp]})\nq = {fmt_p(R.loc[g, 'q'])}, peak {R.loc[g, 'peak_subtype']}", fontsize=6.5)
        ax.tick_params(axis="x", labelsize=6)
        if ax not in axes[:, 0]: ax.set_ylabel("")
    save(fig, "fig8c_boxplots")


if __name__ == "__main__":
    # Figs. 6-8 panels moved to panels_fig6_8.py (30 Oct 2026, Arial 7/6 pt standard). The old fig6*/fig7*/fig8
    # functions below are kept for reference only and are no longer callable from here, so they cannot
    # overwrite the new panels.
    which = sys.argv[2:] or ["5b", "5c", "5c_rho"]
    fn = {"5b": fig5b, "5c": fig5c, "5c_rho": fig5c_rho}
    for w in which: fn[w]()
