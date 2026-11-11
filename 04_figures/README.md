# 04. Figures and supplementary tables

Draws every panel of Figs. 5–8 and Supplementary Fig. 1 as a vector PDF at its exact printed size, lays the
panels out into figure pages, and builds the supplementary tables workbook.

Typography: Arial; 7 pt panel, axis and legend titles; 6 pt tick labels and other small text; 8 pt bold panel letters.

## Run order
| Step | Script | Panels / output | Main inputs |
|---|---|---|---|
| 1 | `fig5_supp_R_panels.R` | Fig. 5a (UMAPs), Fig. 5d (melanoma fraction by segment), Supplementary Fig. 1a (Louvain clusters) and 1b (marker dot plot), marker tables (Supplementary Table 1) | `2026-01-25_melanoma_filter.rds`, GeoMx metadata, `2026-01-28_cibersort_fraction.txt` |
| 2 | `panels.py` | Fig. 5b (signature matrix), Fig. 5c (pseudobulk benchmark) | Signature matrix, pseudobulk `CIBERSORTx_Results.txt` files |
| 3 | `panels_fig6_8.py` | Figs. 6a–e, 7a–f, 8a–c | EcoTyper discovery and recovery outputs, TCGA-SKCM TPM and fractions, GeoMx metadata, HiRes expression, `scan_*.csv` from `03_downstream_statistics/hires_all.py` |
| 4 | `compose.py` | `Figure5.pdf` … `Figure8.pdf`, the combined draft and the supplementary figure page | The panel PDFs from steps 1–3 |
| 5 | `build_supplementary_tables.py` | `Supplementary_Tables.xlsx` (Tables 1–4) | Marker table, `03_downstream_statistics` outputs, EcoTyper recovery output, `scan_*.csv` |

```bash
Rscript fig5_supp_R_panels.R
uv run --with matplotlib --with pandas --with scipy --with numpy python panels.py panels_regenerated
uv run --with matplotlib --with pandas --with scipy --with numpy python panels_fig6_8.py panels_regenerated
python3 compose.py panels_regenerated .
uv run --with pandas --with openpyxl python build_supplementary_tables.py
```

`panels_fig6_8.py` can also draw single panels, for example `python panels_fig6_8.py panels_regenerated 7def`.

## Notes
- **Marker genes:** `fig5_supp_R_panels.R` reuses the saved marker table; set `RERUN_MARKERS=1` to recompute it.
  Markers are from `FindAllMarkers(only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)` on log-normalized
  RNA, kept at adjusted P < 0.05 and ranked by adjusted P, then log2 fold change.
- **Figs. 7d–f** use every melanoma AOI of the compartment, not only the AOIs assigned to an ecotype. Figs. 7a
  and 7b use the 167 assigned AOIs; Fig. 7c labels every melanoma AOI with its higher-abundance ecotype.
- **Fig. 6b** is the only raster element: the 951-gene × 473-tumor heat map is embedded as one image at native
  resolution (one pixel per cell) to keep the PDF small.
- **Seeds:** box-plot point jitter uses a fixed seed (1), so the figures are identical on every run.
- **Software:** R with Seurat 5, dplyr, ggplot2, ggrepel, patchwork, ggtext; Python with matplotlib, pandas,
  numpy, scipy, openpyxl and PyMuPDF (`compose.py`). Arial must be installed.

## Before running
Each script sets the project folder at the top (`base_dir` in R, `BASE` in Python). Change it to your own path.
