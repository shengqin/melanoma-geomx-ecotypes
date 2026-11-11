# Fig. 5 and Supplementary Fig. (scRNA reference) panels drawn in R at their exact printed size.
#
# Typography standard for every panel (author request, 30 Oct 2026):
#   Arial; 7 pt for titles (panel titles, axis titles, legend titles);
#   6 pt for tick labels and all other small text (legend entries, in-plot labels, P values).
# Each PDF is saved at the size it is printed, so compose.py places it unscaled.
#
# Outputs (in Figures/panels_regenerated and Figures/supplementary/source):
#   fig5a_umaps.pdf                 UMAP by study + UMAP by cell type (Fig. 5a)
#   fig5d_melanoma_fraction.pdf     deconvolved melanoma fraction by segment (Fig. 5d)
#   supp_a_umap_clusters.pdf        UMAP by Louvain cluster (Supp. a)
#   SuppFig_scRNA_marker_dotplot.pdf marker-gene dot plot (Supp. b)
#   SuppTable_scRNA_markers_*.csv   marker tables (author marker rule)

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})

# Step 1. Paths and the shared typography.
base_dir <- "/Users/Brian/Library/CloudStorage/Box-Box/NanoString GeoMx data analysis"
figures_dir <- file.path(base_dir, "manuscript/2026-09-15_manuscript/Figures")
panels_dir <- file.path(figures_dir, "panels_regenerated")
# Submission-ready supplementary figures sit in supplementary/; everything they are built from goes in supplementary/source/.
supp_dir <- file.path(figures_dir, "supplementary", "source")
dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)

title_pt <- 7   # titles: panel, axis and legend titles
small_pt <- 6   # tick labels, legend entries, in-plot labels
small_text_mm <- small_pt / ggplot2::.pt   # geom_text sizes are in mm, so convert from pt
pt_to_mm <- 25.4 / 72

# One theme used by every R panel so fonts are identical.
figure_theme <- theme_classic(base_size = title_pt, base_family = "Arial") +
  theme(
    plot.title = element_text(size = title_pt, face = "plain", hjust = 0.5, colour = "black"),
    axis.title = element_text(size = title_pt, colour = "black"),
    axis.text = element_text(size = small_pt, colour = "black"),
    legend.title = element_text(size = title_pt, colour = "black"),
    legend.text = element_text(size = small_pt, colour = "black"),
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    axis.ticks.length = unit(1, "mm"),
    legend.key.size = unit(3, "mm"),
    plot.margin = margin(2, 2, 2, 2, "mm")
  )

save_panel <- function(plot_object, file_path, width_pt, height_pt) {
  # Save at the exact printed size; cairo_pdf embeds Arial.
  ggsave(file_path, plot_object,
         width = width_pt * pt_to_mm, height = height_pt * pt_to_mm,
         units = "mm", device = cairo_pdf)
}

# Step 2. Load the final scRNA reference (the object used to build the signature matrix).
seurat_object <- readRDS(file.path(base_dir, "data/2026-01-25_melanoma_filter.rds"))
DefaultAssay(seurat_object) <- "RNA"

# Cell-type order and display names match Fig. 5b and the Fig. 6a palette.
cell_type_order <- c("Melanoma", "T-CD8", "T-CD4", "NK", "B_cell", "Macrophage", "Endothelial", "CAF", "Mural_cell")
cell_type_names <- c(Melanoma = "Melanoma", `T-CD8` = "T-CD8", `T-CD4` = "T-CD4", NK = "NK", B_cell = "B cell",
                     Macrophage = "Macrophage", Endothelial = "Endothelial", CAF = "CAF", Mural_cell = "Mural cell")
cell_type_colours <- c(Melanoma = "#7B4A23", `T-CD8` = "#2C6FB0", `T-CD4` = "#7FB2E0", NK = "#9C27B0", B_cell = "#3FA34D",
                       Macrophage = "#E08A2B", Endothelial = "#D7263D", CAF = "#C85A54", Mural_cell = "#6D597A")

# Study names and colours are the same as in the original Fig. 5a (ggplot default hues).
study_names <- c(GSE215121_10X = "GSE215121 (10x)", GSE115978_SmartSeq2 = "GSE115978 (Smart-seq2)", `10X_Public` = "Public 10x dataset")
study_colours <- c(GSE215121_10X = "#F8766D", GSE115978_SmartSeq2 = "#00BA38", `10X_Public` = "#619CFF")

# Step 3. Collect UMAP coordinates and labels into one table.
umap_coordinates <- as.data.frame(Embeddings(seurat_object, "umap"))
colnames(umap_coordinates) <- c("UMAP_1", "UMAP_2")
cell_types <- sub("−", "-", as.character(seurat_object$mega_cluster))   # unify minus sign if present
umap_table <- data.frame(
  umap_coordinates,
  study = as.character(seurat_object$Study),
  cell_type = factor(cell_types, levels = cell_type_order),
  cluster = as.character(seurat_object$RNA_snn_res.0.2)
)
stopifnot(!any(is.na(umap_table$cell_type)))

# Draw the large GSE215121 set first so the smaller datasets stay visible on top
# (same drawing order as the original Fig. 5a).
study_draw_order <- c("GSE215121_10X", "GSE115978_SmartSeq2", "10X_Public")
umap_table <- umap_table[order(match(umap_table$study, study_draw_order)), ]

# Points: tiny, borderless dots so 48,733 cells stay readable at print size.
point_size <- 0.05

umap_axes <- list(
  labs(x = "UMAP 1", y = "UMAP 2"),
  coord_equal(),
  figure_theme
)

# Step 4. Fig. 5a, left: UMAP coloured by study, legend below.
study_plot <- ggplot(umap_table, aes(UMAP_1, UMAP_2, colour = study)) +
  geom_point(size = point_size, stroke = 0, shape = 16) +
  scale_colour_manual(values = study_colours, labels = study_names, breaks = study_draw_order, name = NULL) +
  guides(colour = guide_legend(nrow = 2, override.aes = list(size = 1.6))) +
  ggtitle("Study") +
  umap_axes +
  theme(legend.position = "bottom", legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(1, "mm"))

# Step 5. Fig. 5a, right: UMAP coloured by cell type, labelled at each type's median position.
cell_type_label_positions <- umap_table %>%
  group_by(cell_type) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")
cell_type_label_positions$label <- cell_type_names[as.character(cell_type_label_positions$cell_type)]

cell_type_plot <- ggplot(umap_table, aes(UMAP_1, UMAP_2, colour = cell_type)) +
  geom_point(size = point_size, stroke = 0, shape = 16) +
  geom_text_repel(data = cell_type_label_positions, aes(label = label), colour = "black",
                  size = small_text_mm, family = "Arial", seed = 1,          # seed fixes label placement
                  min.segment.length = 0, segment.size = 0.2, box.padding = 0.25) +
  scale_colour_manual(values = cell_type_colours, guide = "none") +
  ggtitle("Cell type") +
  umap_axes

# patchwork aligns the two plotting areas so both UMAPs print at the same size.
fig5a_plot <- study_plot + cell_type_plot + plot_layout(ncol = 2)
# Row 1 of Fig. 5 is 200 pt tall so two square UMAPs plus the study legend fill the 334-pt width.
save_panel(fig5a_plot, file.path(panels_dir, "fig5a_umaps.pdf"), width_pt = 334, height_pt = 200)

# Step 6. Supplementary a: UMAP coloured by Louvain cluster (resolution 0.2), numbers on the plot.
cluster_levels <- as.character(sort(as.integer(unique(umap_table$cluster))))
umap_table$cluster <- factor(umap_table$cluster, levels = cluster_levels)
cluster_label_positions <- umap_table %>%
  group_by(cluster) %>%
  summarise(UMAP_1 = median(UMAP_1), UMAP_2 = median(UMAP_2), .groups = "drop")

cluster_plot <- ggplot(umap_table, aes(UMAP_1, UMAP_2, colour = cluster)) +
  geom_point(size = point_size, stroke = 0, shape = 16) +
  geom_text_repel(data = cluster_label_positions, aes(label = cluster), colour = "black",
                  size = small_text_mm, family = "Arial", seed = 1,
                  min.segment.length = 0, segment.size = 0.2, box.padding = 0.2) +
  scale_colour_discrete(guide = "none") +
  ggtitle("Louvain clusters (resolution 0.2)") +
  umap_axes
save_panel(cluster_plot, file.path(supp_dir, "supp_a_umap_clusters.pdf"), width_pt = 210, height_pt = 200)

# Step 7. Fig. 5d: deconvolved melanoma fraction by segment within each tissue group.
# Same data and test as the original panel (2026-01-15_scRNA.R, "Make fraction plot"):
# Kruskal-Wallis across the three segments, separately for normal skin, nevi and melanoma.
geomx_metadata <- read.table(file.path(base_dir, "data/2025/2025-01-27_GeoMxmetadata.txt"), sep = "\t", header = TRUE)
geomx_fractions <- read.table(file.path(base_dir, "data/2026-01-28_cibersort_fraction.txt"), sep = "\t", header = TRUE)
fraction_table <- merge(geomx_metadata, geomx_fractions, by.x = "ID", by.y = "Mixture")

segment_names <- c(S100 = "S100+", CD45 = "CD45+", `neg/neg` = "S100−/CD45−")
segment_colours <- c(`S100+` = "#619CFF", `CD45+` = "#F8766D", `S100−/CD45−` = "#00BA38")
tissue_names <- c(Norm = "Normal skin", Nevi = "Nevi", Melanoma = "Melanoma")

fraction_table$segment <- factor(segment_names[fraction_table$Segment_Tags], levels = segment_names)
fraction_table$tissue <- factor(tissue_names[fraction_table$Pathophysiology], levels = tissue_names)

kruskal_labels <- data.frame(tissue = factor(tissue_names, levels = tissue_names), label = NA_character_)
for (tissue_index in seq_along(tissue_names)) {
  tissue_rows <- fraction_table[fraction_table$tissue == tissue_names[tissue_index], ]
  kruskal_p <- kruskal.test(Melanoma ~ segment, data = tissue_rows)$p.value
  cat(sprintf("Kruskal-Wallis P, %s: %.3g\n", tissue_names[tissue_index], kruskal_p))
  # Rich-text label (ggtext): everything in Arial, and the exponent kept at 6 pt so no text
  # falls below the 6 pt minimum (Arial has no superscript digits, and plotmath shrinks exponents).
  if (kruskal_p < 0.001) {
    exponent <- floor(log10(kruskal_p))
    mantissa <- sprintf("%.1f", kruskal_p / 10^exponent)
    exponent_text <- sub("-", "\u2212", as.character(exponent))
    kruskal_labels$label[tissue_index] <- paste0("<i>P</i> = ", mantissa, " \u00d7 10<sup><span style='font-size:", small_pt, "pt'>",
                                                 exponent_text, "</span></sup>")
  } else if (kruskal_p < 0.01) {
    kruskal_labels$label[tissue_index] <- paste0("<i>P</i> = ", sprintf("%.3f", kruskal_p))
  } else {
    kruskal_labels$label[tissue_index] <- paste0("<i>P</i> = ", sprintf("%.2f", kruskal_p))
  }
}
top_of_data <- max(fraction_table$Melanoma, na.rm = TRUE)

fraction_plot <- ggplot(fraction_table, aes(tissue, Melanoma, fill = segment)) +
  geom_boxplot(alpha = 0.75, outlier.shape = NA, linewidth = 0.3, width = 0.75,
               position = position_dodge(width = 0.8)) +
  geom_point(position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2, seed = 1),  # seed fixes jitter
             size = 0.4, alpha = 0.7, stroke = 0) +
  ggtext::geom_richtext(data = kruskal_labels, aes(tissue, top_of_data * 1.08, label = label), inherit.aes = FALSE,
                        size = small_text_mm, family = "Arial", fill = NA, label.colour = NA,
                        label.padding = unit(0, "mm")) +
  scale_fill_manual(values = segment_colours, name = "Segment") +
  scale_y_continuous(limits = c(0, top_of_data * 1.14), expand = expansion(mult = c(0.02, 0))) +
  labs(title = "Deconvolved melanoma fraction", x = NULL, y = "Melanoma fraction") +
  figure_theme +
  theme(legend.position = "right")
save_panel(fraction_plot, file.path(panels_dir, "fig5d_melanoma_fraction.pdf"), width_pt = 229, height_pt = 150)

# Step 8. Marker genes (author rule): FindAllMarkers(only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# on log-normalised RNA, no downsampling; rank by p_val_adj then avg_log2FC; keep p_val_adj < 0.05.
# The table from the 21 Oct 2026 run is reused unless RERUN_MARKERS=1.
seurat_object$cell_type <- factor(cell_types, levels = cell_type_order)
Idents(seurat_object) <- "cell_type"
marker_file <- file.path(supp_dir, "SuppTable_scRNA_markers_all.csv")
if (file.exists(marker_file) && Sys.getenv("RERUN_MARKERS") != "1") {
  all_markers <- read.csv(marker_file, check.names = FALSE)
} else {
  if (length(Layers(seurat_object[["RNA"]], search = "counts")) > 1) {
    seurat_object[["RNA"]] <- JoinLayers(seurat_object[["RNA"]])
  }
  all_markers <- FindAllMarkers(seurat_object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
  all_markers <- all_markers %>% arrange(cluster, p_val_adj, desc(avg_log2FC))
  write.csv(all_markers, marker_file, row.names = FALSE)
}

# Report the top-50 QC (how many of the top 50 per type pass p_adj < 0.05).
top50_markers <- all_markers %>% group_by(cluster) %>% slice_head(n = 50) %>% ungroup()
top50_qc <- top50_markers %>% group_by(cluster) %>% summarise(n_top = n(), n_pass = sum(p_val_adj < 0.05), .groups = "drop")
print(top50_qc)
write.csv(top50_markers %>% filter(p_val_adj < 0.05), file.path(supp_dir, "SuppTable_scRNA_markers_top50.csv"), row.names = FALSE)

# Step 9. Choose the dot-plot genes: each significant gene is assigned to the one cell type where it is
# most enriched (so no gene appears twice), then the top 5 per type by p_adj then log2FC.
genes_per_type <- 5
significant_markers <- all_markers %>% filter(p_val_adj < 0.05)
significant_markers$cluster <- factor(as.character(significant_markers$cluster), levels = cell_type_order)
best_type_per_gene <- significant_markers %>% group_by(gene) %>% slice_max(avg_log2FC, n = 1, with_ties = FALSE) %>% ungroup()
chosen_markers <- best_type_per_gene %>%
  arrange(cluster, p_val_adj, desc(avg_log2FC)) %>%
  group_by(cluster) %>%
  slice_head(n = genes_per_type) %>%
  ungroup()
write.csv(chosen_markers, file.path(supp_dir, "SuppFig_dotplot_genes.csv"), row.names = FALSE)
chosen_genes <- chosen_markers$gene

# Step 10. Supplementary b: dot plot. Cell types run top to bottom in Fig. 5 order and genes left to right,
# so each type's genes form a diagonal block.
display_levels <- rev(unname(cell_type_names[cell_type_order]))
seurat_object$cell_type_display <- factor(unname(cell_type_names[cell_types]), levels = display_levels)
dot_plot <- DotPlot(seurat_object, features = chosen_genes, group.by = "cell_type_display",
                    cols = c("grey92", "#B2182B"), dot.scale = 3.6, col.min = -1.5, col.max = 2.5) +
  labs(title = "Marker genes per cell type", x = NULL, y = NULL,
       size = "% cells expressing", colour = "Scaled mean expression") +
  guides(colour = guide_colorbar(order = 1, display = "gradient", barheight = unit(12, "mm"), barwidth = unit(2.5, "mm")),
         size = guide_legend(order = 2)) +
  figure_theme +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, face = "italic", size = small_pt, colour = "black"),
        panel.grid.major = element_line(colour = "grey94", linewidth = 0.25),
        legend.spacing.y = unit(1.5, "mm"))
block_edges <- cumsum(table(factor(chosen_markers$cluster, levels = cell_type_order)))
dot_plot <- dot_plot + geom_vline(xintercept = head(block_edges, -1) + 0.5, colour = "grey70", linewidth = 0.25)
save_panel(dot_plot, file.path(supp_dir, "SuppFig_scRNA_marker_dotplot.pdf"), width_pt = 498, height_pt = 185)

cat("R panels saved\n")
writeLines(capture.output(sessionInfo()), file.path(supp_dir, "sessionInfo_R_panels.txt"))
