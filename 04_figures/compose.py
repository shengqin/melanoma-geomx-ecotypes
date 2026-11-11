"""Compose Figures 5-8 (vector PDF) from regenerated matplotlib panels + R-generated vector panels."""
import sys, os, fitz

BASE = "/Users/Brian/Library/CloudStorage/Box-Box/NanoString GeoMx data analysis"
A = f"{BASE}/analysis"
P = sys.argv[1]            # regenerated panels dir
OUTDIR = sys.argv[2]       # output dir
os.makedirs(OUTDIR, exist_ok=True)

SRC = {
    "umap": f"{A}/2026-01-28_final.pdf",
    "sigmtx": f"{BASE}/manuscript/2026-06-02_manuscript/figure_panels/fig1B_sigmatrix_chessboard.pdf",
    "melfrac": f"{A}/2026-01-28_patho_Melanoma_Fractions.pdf",
    "tcga_eco": f"{A}/2026-01-28_TCGA_SKCM_Discovery_997_4_4/Ecotypes/heatmap_assigned_samples_viridis.pdf",
    "geomx_eco": f"{A}/2026-01-28_GeoMx_recovery_noNorm/Ecotypes/heatmap_assigned_samples_viridis.pdf",
}
# Supplementary source panels (from fig5_supp_R_panels.R) live in supplementary/source/.
SUPP_SOURCE = f"{BASE}/manuscript/2026-09-15_manuscript/Figures/supplementary/source"
SRC["dotplot"] = f"{SUPP_SOURCE}/SuppFig_scRNA_marker_dotplot.pdf"
# Native-size panels from fig5_supp_R_panels.R (drawn at their printed size, Arial 7/6 pt).
SRC["fig5a_umaps"] = f"{P}/fig5a_umaps.pdf"
SRC["fig5d_melanoma_fraction"] = f"{P}/fig5d_melanoma_fraction.pdf"
SRC["supp_a_umap_clusters"] = f"{SUPP_SOURCE}/supp_a_umap_clusters.pdf"

# Native-size Fig. 6-8 panels from panels_fig6_8.py (Arial, 7 pt titles, 6 pt ticks).
for panel_name in ["fig6a_tcga_composition", "fig6b_tcga_melanoma_states", "fig6c_s01_vs_signatures",
                   "fig6d_tcga_ecotype_heatmap", "fig6e_tcga_associations",
                   "fig7a_geomx_ecotype_heatmap", "fig7b_ecotype_by_segment", "fig7c_niche_coupling",
                   "fig7d_S01_by_subtype", "fig7e_stromal_states", "fig7f_immune_states",
                   "fig8a_counts", "fig8b_heatmap", "fig8c_boxplots"]:
    SRC[panel_name] = f"{P}/{panel_name}.pdf"

# Arial font files for page text, so letters and headers match the Arial panels.
ARIAL_FILE = "/System/Library/Fonts/Supplemental/Arial.ttf"
ARIAL_BOLD_FILE = "/System/Library/Fonts/Supplemental/Arial Bold.ttf"
LETTER_PT = 8      # panel letters: 8 pt bold (Nature style)

# The UMAP source PDF (2026-01-28_final.pdf, 1440 x 432 pt) holds three panels side by side.
# These x-ranges were read from the PDF itself (axis-label positions and empty gutters),
# so each panel can be cropped out as vector graphics without re-plotting.
UMAP_STUDY_REGION = (0, 0, 576, 432)        # coloured by study, with its legend
UMAP_CLUSTER_REGION = (576, 0, 1028, 432)   # coloured by Louvain cluster (RNA_snn_res.0.2)
UMAP_CELLTYPE_REGION = (1028, 0, 1440, 432) # coloured by final cell-type annotation

for k in ["fig5b_sigmatrix_full", "fig5c_pseudobulk_scatter", "fig5c_pseudobulk_rho", "fig6a_tcga_composition", "fig6b_tcga_melanoma_states", "fig6c_s01_vs_signatures",
          "fig6e_tcga_associations", "fig7b_S01_by_subtype", "fig7b_ecotype_by_segment", "fig7c_stromal_states", "fig7d_immune_states", "fig7e_niche_coupling",
          "fig8a_counts", "fig8b_heatmap", "fig8c_boxplots"]:
    SRC[k] = f"{P}/{k}.pdf"

PAGE_W, PAGE_H = 612, 792
ML, MT, W = 51, 36, 510      # left margin, top margin (0.5 in; the draft header was removed 30 Oct 2026), content width
GAP, LW = 12, 12          # gap between rows/cols, letter gutter width
_cache = {}


def doc(key):
    if key not in _cache: _cache[key] = fitz.open(SRC[key])
    return _cache[key]


def content_bbox(key, pad=2):
    page = doc(key)[0]; pr = page.rect; bb = None
    def add(r):
        nonlocal bb
        r = fitz.Rect(r)
        if r.is_empty or r.is_infinite: return
        if r.get_area() > 0.9 * pr.get_area(): return        # page-background fill
        bb = r if bb is None else bb | r
    for d in page.get_drawings(): add(d["rect"])
    for b in page.get_text("dict")["blocks"]: add(b["bbox"])
    for img in page.get_images(full=True):
        for r in page.get_image_rects(img[0]): add(r)
    bb = bb & pr
    return fitz.Rect(bb.x0 - pad, bb.y0 - pad, bb.x1 + pad, bb.y1 + pad) & pr


def place(page, key, x, y, width=None, height=None, clip=True, region=None):
    # region: optional (x0, y0, x1, y1) in source-page points; it overrides the automatic crop.
    src = doc(key)
    if region is not None:
        c = fitz.Rect(region)
    elif clip:
        c = content_bbox(key)
    else:
        c = src[0].rect
    ar = c.height / c.width
    if width is None and height is None:
        # Native size: the panel was drawn at its printed size, so place it unscaled.
        width = c.width
        height = c.height
    if width is None: width = height / ar
    if height is None: height = width * ar
    page.show_pdf_page(fitz.Rect(x, y, x + width, y + height), src, 0, clip=c)
    return width, height


def letter(page, s, x, y):
    page.insert_text((x, y + 9), s, fontsize=10, fontname="hebo")


def header(page, s):
    page.insert_text((ML, 26), s, fontsize=8, fontname="helv", color=(0.35, 0.35, 0.35))
    page.insert_text((ML, 36), "Draft assembled 15 Sep 2026 from vector source panels; numbers verified against analysis tables.",
                     fontsize=6.5, fontname="helv", color=(0.55, 0.55, 0.55))


def arial_letter(page, text, x, y):
    # Panel letter in Arial Bold, 8 pt; the baseline sits 7 pt below the panel top.
    page.insert_text((x, y + 7), text, fontsize=LETTER_PT, fontname="ArialBold", fontfile=ARIAL_BOLD_FILE)


def arial_header(page, text):
    # Page title in Arial, matching the panels; same grey draft line as the other pages.
    page.insert_text((ML, 26), text, fontsize=8, fontname="Arial", fontfile=ARIAL_FILE, color=(0.35, 0.35, 0.35))
    page.insert_text((ML, 36), "Draft assembled 15 Sep 2026 from vector source panels; numbers verified against analysis tables.",
                     fontsize=6.5, fontname="Arial", fontfile=ARIAL_FILE, color=(0.55, 0.55, 0.55))


def new_page(out):
    return out.new_page(width=PAGE_W, height=PAGE_H)


# ---------------------------------------------------------------- Figure 5
def figure5(out):
    # Every Fig. 5 panel is drawn at its printed size (Arial, 7 pt titles, 6 pt ticks),
    # so each one is placed unscaled and the fonts stay identical across panels.

    # Step 1. New page (no page header: figure files for submission carry no title text).
    page = new_page(out)

    # Step 2. Row 1: panel a (UMAPs by study and by cell type) and panel b (signature matrix).
    y = MT
    x = ML
    arial_letter(page, "a", x, y)
    x = x + LW
    umap_width, umap_height = place(page, "fig5a_umaps", x, y, clip=False)
    x = x + umap_width + GAP
    arial_letter(page, "b", x, y)
    x = x + LW
    _, heatmap_height = place(page, "fig5b_sigmatrix_full", x, y, clip=False)
    row1_height = max(umap_height, heatmap_height)

    # Step 3. Row 2: panel c (pseudobulk Spearman) and panel d (melanoma fraction by segment).
    y = y + row1_height + GAP
    x = ML
    arial_letter(page, "c", x, y)
    x = x + LW
    rho_width, rho_height = place(page, "fig5c_pseudobulk_rho", x, y, clip=False)
    x = x + rho_width + GAP
    arial_letter(page, "d", x, y)
    x = x + LW
    _, fraction_height = place(page, "fig5d_melanoma_fraction", x, y, clip=False)
    row2_height = max(rho_height, fraction_height)

    # Step 4. Return the bottom of the content so the caller can report page use.
    return y + row2_height


# ---------------------------------------------------------------- Supplementary figure (scRNA reference)
def supplementary_reference_figure(out):
    # Step 1. New page (no page header: figure files for submission carry no title text).
    page = new_page(out)

    # Step 2. Panel a: UMAP coloured by Louvain cluster, placed at its printed size.
    y = MT
    arial_letter(page, "a", ML, y)
    _, cluster_height = place(page, "supp_a_umap_clusters", ML + LW, y, clip=False)

    # Step 3. Panel b: marker-gene dot plot, placed at its printed size.
    y = y + cluster_height + GAP
    arial_letter(page, "b", ML, y)
    _, dotplot_height = place(page, "dotplot", ML + LW, y, clip=False)

    # Step 4. Return the bottom of the content.
    return y + dotplot_height


# ---------------------------------------------------------------- Figures 6-8 (native-size panels)
def place_row(page, y, panels):
    # Place panels left to right at their printed size, each preceded by its letter in a 12-pt gutter.
    # panels: list of (letter, source key). Returns the height of the tallest panel in the row.
    x = ML
    row_height = 0
    for panel_letter, key in panels:
        arial_letter(page, panel_letter, x, y)
        x = x + LW
        panel_width, panel_height = place(page, key, x, y, clip=False)
        x = x + panel_width + GAP
        row_height = max(row_height, panel_height)
    return row_height


def figure6(out):
    # Step 1. New page (no page header: figure files for submission carry no title text).
    page = new_page(out)

    # Step 2. Three rows: a; b + c; d + e.
    y = MT
    y = y + place_row(page, y, [("a", "fig6a_tcga_composition")]) + GAP
    y = y + place_row(page, y, [("b", "fig6b_tcga_melanoma_states"), ("c", "fig6c_s01_vs_signatures")]) + GAP
    y = y + place_row(page, y, [("d", "fig6d_tcga_ecotype_heatmap"), ("e", "fig6e_tcga_associations")])
    return y


def figure7(out):
    # Step 1. New page (no page header: figure files for submission carry no title text).
    page = new_page(out)

    # Step 2. Four rows: a + b; c + d; e; f.
    y = MT
    y = y + place_row(page, y, [("a", "fig7a_geomx_ecotype_heatmap"), ("b", "fig7b_ecotype_by_segment")]) + GAP
    y = y + place_row(page, y, [("c", "fig7c_niche_coupling"), ("d", "fig7d_S01_by_subtype")]) + GAP
    y = y + place_row(page, y, [("e", "fig7e_stromal_states")]) + GAP
    y = y + place_row(page, y, [("f", "fig7f_immune_states")])
    return y


def figure8(out):
    # Step 1. New page (no page header: figure files for submission carry no title text).
    page = new_page(out)

    # Step 2. Two rows: a + b; c.
    y = MT
    y = y + place_row(page, y, [("a", "fig8a_counts"), ("b", "fig8b_heatmap")]) + GAP
    y = y + place_row(page, y, [("c", "fig8c_boxplots")])
    return y


if __name__ == "__main__":
    combined = fitz.open()
    for n, fn in [(5, figure5), (6, figure6), (7, figure7), (8, figure8)]:
        single = fitz.open(); ymax = fn(single)
        print(f"Figure {n}: content bottom y={ymax:.0f} of {PAGE_H}")
        single.save(f"{OUTDIR}/Figure{n}.pdf", garbage=3, deflate=True)
        combined.insert_pdf(single)
        single[0].get_pixmap(dpi=110).save(f"{OUTDIR}/Figure{n}_preview.png")
    combined.save(f"{OUTDIR}/Figures5-8_draft.pdf", garbage=3, deflate=True)

    # The supplementary figure for submission goes in supplementary/; its preview goes in supplementary/source/.
    supp_dir = os.path.join(OUTDIR, "supplementary")
    supp_source_dir = os.path.join(supp_dir, "source")
    os.makedirs(supp_source_dir, exist_ok=True)
    supp = fitz.open()
    supp_bottom = supplementary_reference_figure(supp)
    print(f"Supplementary reference figure: content bottom y={supp_bottom:.0f} of {PAGE_H}")
    supp.save(os.path.join(supp_dir, "SuppFig_scRNA_reference.pdf"), garbage=3, deflate=True)
    supp[0].get_pixmap(dpi=110).save(os.path.join(supp_source_dir, "SuppFig_scRNA_reference_preview.png"))
    print("done")
