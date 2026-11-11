library(Seurat)
library(dplyr)
library(harmony)
library(ggplot)
library(ggpubr)
file_path <- "/Users/Brian/Library/CloudStorage/Box-Box/NanoString\ GeoMx\ data\ analysis/data"
setwd(file_path)
meta_atlas = readRDS(paste0(file_path,"/2025-10-14_scRNA_melanoma.rds"))
meta_atlas <- FindClusters(meta_atlas, resolution = 0.2) # Start with 0.5, adjust if needed


p_vin = VlnPlot(meta_atlas, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

meta_atlas$cell_type_annotation = case_when(
  meta_atlas$RNA_snn_res.0.05 %in% c(2,3,4,5,6,7,9,10,13,14,15,16) ~ "Melanoma",
  meta_atlas$RNA_snn_res.0.05 == 0 ~ "T_NK",
  meta_atlas$RNA_snn_res.0.05 == 1 ~ "B_cell",
  meta_atlas$RNA_snn_res.0.05 == 8 ~ "CAF",
  meta_atlas$RNA_snn_res.0.05 == 11 ~ "Endothelial",
  meta_atlas$RNA_snn_res.0.05 == 12 ~ "Macrophage",
)

p1 <- DimPlot(meta_atlas, reduction = "umap", group.by = c("Study", "cell_type_annotation", "RNA_snn_res.0.05"), raster = FALSE, label = T) + 
  ggtitle("Integration by Study")

## Subcluster for melanoma cell
meta_melanoma = subset(meta_atlas, cell_type_annotation == "Melanoma")
meta_melanoma = recluster(meta_melanoma)
meta_melanoma = FindClusters(meta_melanoma, resolution = 0.2)

meta_melanoma_markers = find_target_markers(meta_melanoma, c(15,18,21,22,23,24))

meta_melanoma2 = subset(meta_melanoma, RNA_snn_res.0.2 %in% c(15,23,24), invert = T)
meta_melanoma2 = recluster(meta_melanoma2)
meta_melanoma2 = FindClusters(meta_melanoma2, resolution = 0.2)
p_melanoma2 <- DimPlot(meta_melanoma2, reduction = "umap", group.by = c("Study", "RNA_snn_res.0.2"), raster = FALSE, label = T) + 
  ggtitle("Integration by Study")

meta_melanoma3 = subset(meta_melanoma2, RNA_snn_res.0.2 %in% c(22,24), invert = T)
meta_melanoma3 = recluster(meta_melanoma3)
meta_melanoma3 = FindClusters(meta_melanoma3, resolution = 0.2)


meta_melanoma4 = subset(meta_melanoma3, RNA_snn_res.0.2 %in% c(9), invert = T)
meta_melanoma4 = recluster(meta_melanoma4)
meta_melanoma4 = FindClusters(meta_melanoma4, resolution = 0.2)


meta_melanoma5 = subset(meta_melanoma4, RNA_snn_res.0.2 %in% c(6,16,19), invert = T)
meta_melanoma5 = recluster(meta_melanoma5)
meta_melanoma5 = FindClusters(meta_melanoma5, resolution = 0.2)


meta_melanoma6 = subset(meta_melanoma5, RNA_snn_res.0.2 %in% c(10,17,18,19), invert = T)
meta_melanoma6 = recluster(meta_melanoma6)
meta_melanoma6 = FindClusters(meta_melanoma6, resolution = 0.2)

meta_melanoma7 = subset(meta_melanoma6, RNA_snn_res.0.2 %in% c(15,16), invert = T)
meta_melanoma7 = recluster(meta_melanoma7)
meta_melanoma7 = FindClusters(meta_melanoma7, resolution = 0.2)

meta_melanoma8 = subset(meta_melanoma7, RNA_snn_res.0.2 %in% c(11), invert = T)
meta_melanoma8 = recluster(meta_melanoma8)
meta_melanoma8 = FindClusters(meta_melanoma8, resolution = 0.2)

p_melanoma8 <- DimPlot(
  meta_melanoma8, 
  reduction = "umap", 
  group.by = c("Study", "RNA_snn_res.0.2"), 
  
  # This puts 10X_Public and SmartSeq2 ON TOP.
  # GSE215121_10X will be at the bottom by default.
  order = c("10X_Public", "GSE115978_SmartSeq2"), 
  
  raster = FALSE, 
  label = TRUE
) + ggtitle("Integration by Study")
meta_melanoma_markers = find_target_markers(meta_melanoma8, c(13))
write.csv(meta_melanoma_markers, paste0(file_path,"/melanoma_clusters.csv"), row.names = FALSE)


saveRDS(meta_melanoma8, paste0(file_path,"/2026-01-15_melanoma_cell.rds"))


### Annotate T-NK cluster
meta_tnk = subset(meta_atlas, cell_type_annotation == "T_NK")
meta_tnk = recluster(meta_tnk)
meta_tnk = FindClusters(meta_tnk, resolution = 0.2)

meta_tnk2 = subset(meta_tnk, RNA_snn_res.0.2 %in% c(4,9), invert = T)
meta_tnk2 = recluster(meta_tnk2)
meta_tnk2 = FindClusters(meta_tnk2, resolution = 0.2)

meta_tnk3 = subset(meta_tnk2, RNA_snn_res.0.2 %in% c(6), invert = T)
meta_tnk3 = recluster(meta_tnk3)
meta_tnk3 = FindClusters(meta_tnk3, resolution = 0.2)

meta_tnk4 = subset(meta_tnk3, RNA_snn_res.0.2 %in% c(4), invert = T)
meta_tnk4 = recluster(meta_tnk4)
meta_tnk4 = FindClusters(meta_tnk4, resolution = 0.2)

p_meta_tnk <- DimPlot(
  meta_tnk4, 
  reduction = "umap", 
  group.by = c("Study", "RNA_snn_res.0.2"), 
  
  # This puts 10X_Public and SmartSeq2 ON TOP.
  # GSE215121_10X will be at the bottom by default.
  order = c("10X_Public", "GSE115978_SmartSeq2"), 
  
  raster = FALSE, 
  label = TRUE
) + ggtitle("Integration by Study") + NoLegend()

meta_tnk_markers = find_markers(meta_tnk4)
write.csv(meta_tnk_markers, paste0(file_path,"/melanoma_clusters.csv"), row.names = FALSE)

meta_tnk4$cell_type_annotation = case_when(
  meta_tnk4$RNA_snn_res.0.2 %in% c(0) ~ "T-Naive",
  meta_tnk4$RNA_snn_res.0.2 %in% c(2) ~ "T-CD8",
  meta_tnk4$RNA_snn_res.0.2 == 1 ~ "NK",
  meta_tnk4$RNA_snn_res.0.2 == 4 ~ "T-ISG",
  meta_tnk4$RNA_snn_res.0.2 %in% c(3,5) ~ "Melanoma"
)



saveRDS(meta_tnk4, paste0(file_path,"/2026-01-15_TNK_cell.rds"))

### Annotate B-Macro cluster
meta_bmacro = subset(meta_atlas7, RNA_snn_res.0.2 %in% c(13))
meta_bmacro = recluster(meta_bmacro)
meta_bmacro = FindClusters(meta_bmacro, resolution = 0.2)

meta_bmacro$CD4_score = FetchData(meta_bmacro, vars = "CD4", layer = "data")$CD4
meta_bmacro$CD68_score = FetchData(meta_bmacro, vars = "CD68", layer = "data")$CD68
meta_bmacro$CD468 = meta_bmacro$CD4_score-meta_bmacro$CD68_score
meta_bmacro2 = subset(meta_bmacro, RNA_snn_res.0.2 %in% c(5), invert = T)
meta_bmacro2 = recluster(meta_bmacro2)
meta_bmacro2 = FindClusters(meta_bmacro2, resolution = 0.2)

meta_bmacro3 = subset(meta_bmacro2, RNA_snn_res.0.2 %in% c(7), invert = T)
meta_bmacro3 = recluster(meta_bmacro3)
meta_bmacro3 = FindClusters(meta_bmacro3, resolution = 0.2)

meta_bmacro4 = subset(meta_bmacro3, RNA_snn_res.0.2 %in% c(7), invert = T)
meta_bmacro4 = recluster(meta_bmacro4)
meta_bmacro4 = FindClusters(meta_bmacro4, resolution = 0.2)

p_meta_bmacro <- DimPlot(
  meta_bmacro2, 
  reduction = "umap", 
  group.by = c("mega_cluster", "RNA_snn_res.0.2"), 
  
  # This puts 10X_Public and SmartSeq2 ON TOP.
  # GSE215121_10X will be at the bottom by default.
  order = c("10X_Public", "GSE115978_SmartSeq2"), 
  
  raster = FALSE, 
  label = TRUE
) + ggtitle("Integration by Study") + NoLegend()

meta_bmacro_markers = find_markers(meta_bmacro4)
meta_bmacro_markers = find_target_markers(meta_bmacro3,c(7))
write.csv(meta_bmacro_markers, paste0(file_path,"/melanoma_clusters.csv"), row.names = FALSE)
meta_bmacro4$cell_type_annotation = case_when(
  meta_bmacro4$RNA_snn_res.0.1 %in% c(0) ~ "B-cell",
  meta_bmacro4$RNA_snn_res.0.1 %in% c(1) ~ "Macrophage",
  meta_bmacro4$RNA_snn_res.0.1 == 2 ~ "T-CD8",
  meta_bmacro4$RNA_snn_res.0.1 == 3 ~ "pDC",
  meta_bmacro4$RNA_snn_res.0.1 %in% c(4) ~ "Plasma"
)
saveRDS(meta_bmacro4, paste0(file_path,"/2026-01-15_B-Macro.rds"))

### Annotate fibro-endo
meta_fibend = subset(meta_atlas, cell_type_annotation %in% c("Endothelial","CAF"))
meta_fibend = recluster(meta_fibend)
meta_fibend = FindClusters(meta_fibend, resolution = 0.1)

meta_fibend2 = subset(meta_fibend, RNA_snn_res.0.1 %in% c(7,8,10,11), invert = T)
meta_fibend2 = recluster(meta_fibend2)
meta_fibend2 = FindClusters(meta_fibend2, resolution = 0.2)

meta_fibend3 = subset(meta_fibend2, RNA_snn_res.0.2 %in% c(2,10,12,9), invert = T)
meta_fibend3 = recluster(meta_fibend3)
meta_fibend3 = FindClusters(meta_fibend3, resolution = 0.2)

meta_fibend4 = subset(meta_fibend3, RNA_snn_res.0.2 %in% c(11), invert = T)
meta_fibend4 = recluster(meta_fibend4)
meta_fibend4 = FindClusters(meta_fibend4, resolution = 0.2)

meta_fibend5 = subset(meta_fibend4, RNA_snn_res.0.2 %in% c(11,13), invert = T)
meta_fibend5 = recluster(meta_fibend5)
meta_fibend5 = FindClusters(meta_fibend5, resolution = 0.2)

p_meta_fibend <- DimPlot(
  meta_fibend5, 
  reduction = "umap", 
  group.by = c("Study", "RNA_snn_res.0.2"), 
  
  # This puts 10X_Public and SmartSeq2 ON TOP.
  # GSE215121_10X will be at the bottom by default.
  order = c("10X_Public", "GSE115978_SmartSeq2"), 
  
  raster = FALSE, 
  label = TRUE
) + ggtitle("Integration by Study") + NoLegend()

meta_fibend_markers = find_target_markers(meta_fibend5, c(8,9,10,11))
meta_fibend_markers = find_markers(meta_fibend5)


write.csv(meta_fibend_markers, paste0(file_path,"/melanoma_clusters.csv"), row.names = FALSE)
meta_fibend5$cell_type_annotation = case_when(
  meta_fibend5$RNA_snn_res.0.2 %in% c(2,4,8,7,9) ~ "CAF",
  meta_fibend5$RNA_snn_res.0.2 %in% c(1,3,5) ~ "Endothelial",
  meta_fibend5$RNA_snn_res.0.2 %in% c(0,6) ~ "Mural-cell",
)
saveRDS(meta_fibend5, paste0(file_path,"/2026-01-15_fibro-endo.rds"))

#### Filter

all_clusters <- c(meta_fibend5$cell_type_annotation,
                  meta_melanoma8$cell_type_annotation,
                  meta_bmacro4$cell_type_annotation,
                  meta_tnk4$cell_type_annotation)

# 2. Create the dataframe (Row names are automatically preserved from the vector names)
final_metadata <- data.frame(mega_cluster = all_clusters)
meta_atlas2 = subset(meta_atlas, cells = names(all_clusters))
meta_atlas2 <- AddMetaData(object = meta_atlas2, metadata = final_metadata)
meta_atlas2 = recluster(meta_atlas2)
meta_atlas2 <- FindClusters(meta_atlas2, resolution = 0.2)
p_meta_atlas2 = DimPlot(meta_atlas2, group.by = c("Study","RNA_snn_res.0.2", "mega_cluster"),
                        order = c("10X_Public", "GSE115978_SmartSeq2"),
                        reduction = "umap", label = TRUE) + NoLegend()

saveRDS(meta_atlas2, paste0(file_path,"/2026-01-15_melanoma_scRNA_filter.rds"))

meta_atlas2 = readRDS(paste0(file_path,"/2026-01-15_melanoma_scRNA_filter.rds"))
meta_atlas2$major_cluster1 = case_when(
  meta_atlas2$RNA_snn_res.0.2 %in% c(11,12,9) ~ "caf-endo",
  meta_atlas2$RNA_snn_res.0.2 %in% c(2,5,4,3,8,7,6) ~ "melanoma cell",
  meta_atlas2$RNA_snn_res.0.2 %in% c(1,10,13) ~ "b-macro",
  meta_atlas2$RNA_snn_res.0.2 %in% c(0) ~ "t-nk"
)

meta_atlas2$major_cluster2 <- case_when(
  meta_atlas2$mega_cluster %in% c("B-cell", "Plasma", "pDC", "Macrophage") ~ "b-macro",
  meta_atlas2$mega_cluster %in% c("T-Naive", "T-CD8", "T-ISG", "NK") ~ "t-nk",
  meta_atlas2$mega_cluster %in% c("Melanoma") ~ "melanoma cell",
  meta_atlas2$mega_cluster %in% c("CAF", "Mural-cell", "Endothelial") ~ "caf-endo"
)

DimPlot(meta_atlas2, group.by = c("major_cluster1","major_cluster2"), reduction = "umap", label = TRUE)

meta_atlas2@meta.data$mega_merge = NA
meta_atlas2@meta.data$mega_merge <- apply(meta_atlas2@meta.data %>% select(major_cluster2, major_cluster1), 1, function(x) {
  if (max(table(x)) > 1) {
    # Get the most common annotation
    most_common <- names(which.max(table(x)))
    return(most_common)
  } else {
    return("not defined")
  }
})

meta_atlas3 = subset(meta_atlas2, mega_merge == "not defined", invert = T)
meta_atlas3 = recluster(meta_atlas3)

meta_atlas3 <- FindClusters(meta_atlas3, resolution = 0.2)
p_meta_atlas3 = DimPlot(meta_atlas3, group.by = c("Study","RNA_snn_res.0.2", "mega_cluster"),
                        order = c("10X_Public", "GSE115978_SmartSeq2"),
                        reduction = "umap", label = TRUE) + NoLegend()

meta_atlas3$major_cluster1 = case_when(
  meta_atlas3$RNA_snn_res.0.2 %in% c(8,12,9) ~ "caf-endo",
  meta_atlas3$RNA_snn_res.0.2 %in% c(3,7,13,11,6,5,4,2,16) ~ "melanoma cell",
  meta_atlas3$RNA_snn_res.0.2 %in% c(1,10,15) ~ "b-macro",
  meta_atlas3$RNA_snn_res.0.2 %in% c(0,14) ~ "t-nk"
)

meta_atlas3@meta.data$mega_merge = NA
meta_atlas3@meta.data$mega_merge <- apply(meta_atlas3@meta.data %>% select(major_cluster2, major_cluster1), 1, function(x) {
  if (max(table(x)) > 1) {
    # Get the most common annotation
    most_common <- names(which.max(table(x)))
    return(most_common)
  } else {
    return("not defined")
  }
})
DimPlot(meta_atlas3, group.by = "mega_merge")

saveRDS(meta_atlas3, paste0(file_path,"/2026-01-15_melanoma_scRNA_filter.rds"))

meta_atlas4 = subset(meta_atlas3, mega_merge == "not defined", invert = T)
meta_atlas4 = recluster(meta_atlas4)

meta_atlas4 <- FindClusters(meta_atlas4, resolution = 0.2)
p_meta_atlas4 = DimPlot(meta_atlas4, group.by = c("Study","RNA_snn_res.0.2", "mega_cluster"),
                        order = c("10X_Public", "GSE115978_SmartSeq2"),
                        reduction = "umap", label = TRUE) + NoLegend()

meta_atlas4$major_cluster1 = case_when(
  meta_atlas4$RNA_snn_res.0.2 %in% c(8,12,9) ~ "caf-endo",
  meta_atlas4$RNA_snn_res.0.2 %in% c(3,7,14,11,6,5,4,2,16) ~ "melanoma cell",
  meta_atlas4$RNA_snn_res.0.2 %in% c(1,10,15) ~ "b-macro",
  meta_atlas4$RNA_snn_res.0.2 %in% c(0,13) ~ "t-nk"
)

meta_atlas4@meta.data$mega_merge = NA
meta_atlas4@meta.data$mega_merge <- apply(meta_atlas4@meta.data %>% select(major_cluster2, major_cluster1), 1, function(x) {
  if (max(table(x)) > 1) {
    # Get the most common annotation
    most_common <- names(which.max(table(x)))
    return(most_common)
  } else {
    return("not defined")
  }
})
DimPlot(meta_atlas4, group.by = "mega_merge")

meta_atlas5 <- subset(meta_atlas4, subset = mega_merge != "not defined" & RNA_snn_res.0.2 != "16")
meta_atlas5 = recluster(meta_atlas5)

meta_atlas5 <- FindClusters(meta_atlas5, resolution = 0.2)
p_meta_atlas5 = DimPlot(meta_atlas5, group.by = c("Study","RNA_snn_res.0.2", "mega_cluster"),
                        order = c("10X_Public", "GSE115978_SmartSeq2"),
                        reduction = "umap", label = TRUE) + NoLegend()

saveRDS(meta_atlas5, paste0(file_path,"/2026-01-15_melanoma_scRNA_filter.rds"))

meta_16_markers = find_target_markers(meta_atlas5, c(16))

### Reannotate T-NK
meta_tnk = subset(meta_atlas5, RNA_snn_res.0.2 %in% c(0,12))
meta_tnk = recluster(meta_tnk)
meta_tnk = FindClusters(meta_tnk, resolution = 0.6)

meta_tnk2 = subset(meta_tnk, subset = IGK_Score1 < 0)
meta_tnk2 = recluster(meta_tnk2)
meta_tnk2 = FindClusters(meta_tnk2, resolution = 1.4)

meta_tnk3 = subset(meta_tnk2, subset = RNA_snn_res.1.4 != 2)
meta_tnk3 = recluster(meta_tnk3)
meta_tnk3 = FindClusters(meta_tnk3, resolution = 1.4)

meta_tnk4 = subset(meta_tnk3, subset = RNA_snn_res.1.4 != 14)
meta_tnk4 = recluster(meta_tnk4)
meta_tnk4 = FindClusters(meta_tnk4, resolution = 1.8)

meta_tnk5 = subset(meta_tnk4, subset = RNA_snn_res.1.8 != 9)
meta_tnk5 = recluster(meta_tnk5)
meta_tnk5 = FindClusters(meta_tnk5, resolution = 1.8)

meta_tnk6 = subset(meta_tnk5, subset = RNA_snn_res.1.8 != 0)
meta_tnk6 = recluster(meta_tnk6)
meta_tnk6 = FindClusters(meta_tnk6, resolution = 1.6)

meta_tnk7 = subset(meta_tnk6, subset = RNA_snn_res.1.6 != 14)
meta_tnk7 = recluster(meta_tnk7)
meta_tnk7 = FindClusters(meta_tnk7, resolution = 1.6)

meta_tnk8 = subset(meta_tnk7, subset = RNA_snn_res.1.6 != 0)
meta_tnk8 = recluster(meta_tnk8)
meta_tnk8 = FindClusters(meta_tnk8, resolution = 1.6)

meta_tnk9 = subset(meta_tnk8, subset = RNA_snn_res.1.6 != 15)
meta_tnk9 = recluster(meta_tnk9)
meta_tnk9 = FindClusters(meta_tnk9, resolution = 1.6)

meta_tnk10 = subset(meta_tnk9, subset = percent.trash_loc < 1)
meta_tnk10 = recluster(meta_tnk10)
meta_tnk10 = FindClusters(meta_tnk10, resolution = 1.6)

meta_tnk11 = subset(meta_tnk10, subset = percent.stress < 25)
meta_tnk11 = recluster(meta_tnk11)
meta_tnk11 = FindClusters(meta_tnk11, resolution = 0.1)

meta_tnk12 = subset(meta_tnk11, subset = RNA_snn_res.0.2 == 0 & umap_1 > -1, invert = T)
meta_tnk12 = recluster(meta_tnk12)
meta_tnk12 = FindClusters(meta_tnk12, resolution = 2)


meta_tnk12$mega_cluster = case_when(
  meta_tnk12$RNA_snn_res.2 %in% c(7,5,9,2,12,0,6) ~ "T-CD4",
  meta_tnk12$RNA_snn_res.2 %in% c(14,11,8,13,3) ~ "T-CD8",
  meta_tnk12$RNA_snn_res.2 %in% c(1,10,15,4,16) ~ "NK"
)

cell_remove = setdiff(colnames(meta_tnk),colnames(meta_tnk12))
meta_atlas6 <- subset(meta_atlas5, cells = cell_remove, invert = TRUE)
meta_atlas6@meta.data[colnames(meta_tnk12), "mega_cluster"] <- as.character(meta_tnk12$mega_cluster)
meta_atlas6 = recluster(meta_atlas6)
meta_atlas6 = FindClusters(meta_atlas6, resolution = 0.2)
saveRDS(meta_atlas6, paste0(file_path,"/2026-01-17_melanoma_scRNA_filter.rds"))
p_meta = DimPlot(meta_atlas6, group.by = c("RNA_snn_res.0.2", "mega_cluster"),
                     order = c("10X_Public", "GSE115978_SmartSeq2"),
                     reduction = "umap", label = TRUE) 
meta_tnk_markers = find_target_markers(meta_tnk12, c(14))

trash_genes <- c(
  "LOC646214", "LOC643406", "LOC100131257", "LOC100190986", "LOC286437",
  "LOC100506190", "LOC284454", "LOC284023", "LOC100289019", "LOC100130231",
  "LOC100506385", "LOC100128288", "LOC100128531", "LOC100130581",
  "LOC100505876", "LOC729987", "LOC100506472"
)
stress_genes <- c(
  "MALAT1", "NEAT1",          # Nuclear lncRNAs (high in stripped nuclei)
  "MT1X", "MT1E", "MT2A",     # Metallothioneins (dissociation artifact)
  "HSPA6",                    # Extreme heat shock (terminal stress)
  "LOC100131257", "LOC643406" # Your specific 'trash' LOC genes
)
meta_tnk9[["percent.stress"]] <- PercentageFeatureSet(meta_tnk9, features = stress_genes)

meta_tnk9[["percent.trash_loc"]] <- PercentageFeatureSet(meta_tnk9, features = trash_genes)


meta_tnk_markers = find_target_markers(meta_tnk12, 4)
meta_tnk_markers = find_markers(meta_tnk12)
p_vin = VlnPlot(meta_tnk2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

igk_genes <- grep(pattern = "^IGK", x = rownames(meta_tnk2), value = TRUE)
meta_tnk <- AddModuleScore(
  object = meta_tnk,
  features = list(igk_genes),
  name = "IGK_Score" 
)

p_vin = VlnPlot(meta_tnk, features = "IGK_Score1", pt.size = 0) + 
  geom_boxplot(width=0.1, fill="white", outlier.shape = NA) +
  ggtitle("IGK Contamination Score by Cluster")

write.csv(meta_tnk_markers, paste0(file_path,"/melanoma_clusters.csv"), row.names = FALSE)

meta_atlas7 = readRDS(paste0(file_path,"/2026-01-17_melanoma_filter.rds"))
# 1. Define your markers in a specific order matching the cell types
#    This dictionary structure helps you keep track of which gene is for which cell
marker_list <- list(
  "Melanoma"    = c("PMEL", "MLANA", "MITF"),
  "T-CD8"       = c("CD8A", "CD8B", "GZMB"),
  "T-CD4"       = c("CD3D", "IL7R", "FOXP3","CCR7"),
  "NK"          = c("GNLY", "NKG7", "KLRB1"),
  "B-cell"      = c("MS4A1", "CD79A", "CD19"),
  "Macrophage"  = c("CD68", "CD163", "AIF1"),
  "Endothelial" = c("VWF", "PECAM1", "PLVAP"),
  "Mural−cell"  = c("RGS5", "ACTA2", "MCAM"),
  "CAF"         = c("COL1A1", "DCN", "FAP")
)

# Flatten the list into a simple vector of genes for Seurat
genes_to_plot <- unique(unlist(marker_list))

# 2. Set the order of your clusters (Y-axis)
#    This ensures they appear in the same order as your gene list above
desired_order <- c("Melanoma", "T-CD8", "T-CD4", "NK", "B-cell", 
                   "Macrophage", "Endothelial", "Mural−cell", "CAF")

# Apply this order to the Seurat object
meta_atlas7$mega_cluster <- factor(x = meta_atlas7$mega_cluster, levels = desired_order)

# 3. Generate the organized DotPlot
DotPlot(meta_atlas7, features = genes_to_plot, group.by = "mega_cluster") + 
  RotatedAxis() +                       # Rotate x-axis labels
  coord_flip() +                        # OPTIONAL: Flip to put Cell Types on X-axis, Genes on Y
  #scale_color_distiller(palette = "RdBu") +   # Nice Red-Blue color scheme
  ggtitle("Cell Type Markers") +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))
saveRDS(meta_tnk12, file = paste0(file_path,"/2026-01-17_TNK_cell.rds"))

### 3rd time clean

meta_atlas6$major_cluster1 = case_when(
  meta_atlas6$RNA_snn_res.0.2 %in% c(8,12,9) ~ "caf-endo",
  meta_atlas6$RNA_snn_res.0.2 %in% c(6,11,3,7,13,4,2,5) ~ "melanoma cell",
  meta_atlas6$RNA_snn_res.0.2 %in% c(1,10,14) ~ "b-macro",
  meta_atlas6$RNA_snn_res.0.2 %in% c(0) ~ "t-nk"
)

meta_atlas6$major_cluster2 <- case_when(
  meta_atlas6$mega_cluster %in% c("B-cell", "Plasma", "pDC", "Macrophage") ~ "b-macro",
  meta_atlas6$mega_cluster %in% c("T-CD4", "T-CD8", "NK") ~ "t-nk",
  meta_atlas6$mega_cluster %in% c("Melanoma") ~ "melanoma cell",
  meta_atlas6$mega_cluster %in% c("CAF", "Mural-cell", "Endothelial") ~ "caf-endo",
  meta_atlas6$mega_cluster %in% c("T-Naive", "T-ISG") ~ "not defined",
)

DimPlot(meta_atlas6, group.by = c("major_cluster1","major_cluster2"), reduction = "umap", label = TRUE)

meta_atlas6@meta.data$mega_merge = NA
meta_atlas6@meta.data$mega_merge <- apply(meta_atlas6@meta.data %>% select(major_cluster2, major_cluster1), 1, function(x) {
  if (max(table(x)) > 1) {
    # Get the most common annotation
    most_common <- names(which.max(table(x)))
    return(most_common)
  } else {
    return("not defined")
  }
})

b_macro_diff = setdiff(colnames(meta_bmacro),colnames(meta_bmacro2))
meta_atlas7 = subset(meta_atlas7, mega_cluster %in% c("Plasma", "pDC"), invert = T)
meta_atlas7 = recluster(meta_atlas7)
meta_atlas7 = FindClusters(meta_atlas7, resolution = 0.2)

meta_atlas7$mega_cluster[meta_atlas7$RNA_snn_res.0.2 == "10"] <- "CAF"
meta_atlas7$mega_cluster[meta_atlas7$RNA_snn_res.0.2 == "12"] <- "Mural−cell"
meta_atlas7$mega_cluster[meta_atlas7$RNA_snn_res.0.2 == "9"] <- "Endothelial"
meta_atlas7$mega_cluster[meta_atlas7$RNA_snn_res.0.2 == "13"] <- "Macrophage"
meta_atlas7$mega_cluster[meta_atlas7$RNA_snn_res.0.2 == "2"] <- "B-cell"
meta_atlas7$mega_cluster[meta_atlas7$RNA_snn_res.0.2 == "3"] <- "Melanoma"


remove_cells <- WhichCells(meta_atlas7, expression = RNA_snn_res.0.2 == "3" & mega_cluster == "T-CD4")
meta_atlas7$mega_cluster[pdc_cells] <- "pDC"

macro_cells <- WhichCells(meta_bmacro2, expression = RNA_snn_res.0.2 %in% c(3,4))
meta_atlas7$mega_cluster[macro_cells] <- "Macrophage"

b_cells <- WhichCells(meta_bmacro2, expression = RNA_snn_res.0.2 %in% c(0,1,2))
meta_atlas7$mega_cluster[b_cells] <- "B-cell"

p_meta = DimPlot(meta_atlas7, group.by = c("Study","RNA_snn_res.0.2", "mega_cluster"),
                 order = c("10X_Public", "GSE115978_SmartSeq2"),
                 reduction = "umap", label = TRUE) 

meta_atlas7$major_cluster1 = case_when(
  meta_atlas7$RNA_snn_res.0.2 %in% c(8,12,9) ~ "caf-endo",
  meta_atlas7$RNA_snn_res.0.2 %in% c(6,11,3,7,13,4,2,5) ~ "melanoma cell",
  meta_atlas7$RNA_snn_res.0.2 %in% c(1,10,14) ~ "b-macro",
  meta_atlas7$RNA_snn_res.0.2 %in% c(0) ~ "t-nk"
)

meta_atlas7@meta.data$mega_merge = NA
meta_atlas7@meta.data$mega_merge <- apply(meta_atlas7@meta.data %>% select(major_cluster2, major_cluster1), 1, function(x) {
  if (max(table(x)) > 1) {
    # Get the most common annotation
    most_common <- names(which.max(table(x)))
    return(most_common)
  } else {
    return("not defined")
  }
})
DimPlot(meta_atlas7, group.by = c("mega_merge"), reduction = "umap", label = TRUE)

meta_markers = find_target_markers(meta_atlas7, 12)
write.csv(meta_fibend_markers, paste0(file_path,"/melanoma_clusters.csv"), row.names = FALSE)
saveRDS(meta_atlas7, file = paste0(file_path,"/2026-01-17_melanoma_filter.rds"))

### Refine-melanoma cluster ----
meta_atlas7 = readRDS("../data/2026-01-17_melanoma_filter.rds")

meta_melanoma_ori = subset(meta_atlas7, RNA_snn_res.0.2 %in% c(6,10,8,7,1,3,4,5))
meta_melanoma = recluster(meta_melanoma)
meta_melanoma = FindClusters(meta_melanoma, resolution = 0.2)

meta_melanoma = subset(meta_melanoma, TNK_SigScore1 < 0.02)
meta_melanoma = subset(meta_melanoma, Trash_Score1 < 0.02)
meta_melanoma = subset(meta_melanoma, T_Score1 < 0.05)
meta_melanoma = subset(meta_melanoma, RNA_snn_res.0.2 != 11)
meta_melanoma = recluster(meta_melanoma)
meta_melanoma = FindClusters(meta_melanoma, resolution = 0.2)

b_mela_diff = setdiff(colnames(meta_melanoma_ori),colnames(meta_melanoma))
meta_atlas8 = subset(meta_atlas7, cells = b_mela_diff, invert = T)
meta_atlas8 = recluster(meta_atlas8)
meta_atlas8 = FindClusters(meta_atlas8, resolution = 0.2)
meta_atlas9 = subset(meta_atlas8, mega_merge != "not defined")
meta_atlas9 = recluster(meta_atlas9)
meta_atlas9 = FindClusters(meta_atlas9, resolution = 0.2)

meta_atlas9$mega_cluster[meta_atlas9$RNA_snn_res.0.2 == "10"] <- "CAF"
meta_atlas9$mega_cluster[meta_atlas9$RNA_snn_res.0.2 == "11"] <- "Mural_cell"
meta_atlas9$mega_cluster[meta_atlas9$RNA_snn_res.0.2 == "9"] <- "Endothelial"
meta_atlas9$mega_cluster[meta_atlas9$RNA_snn_res.0.2 == "12"] <- "Macrophage"
meta_atlas9$mega_cluster[meta_atlas9$RNA_snn_res.0.2 == "1"] <- "B_cell"
meta_atlas9$mega_cluster[meta_atlas9$RNA_snn_res.0.2 %in% c(5,4,0,2,3,8,13)] <- "Melanoma"
saveRDS(meta_atlas9, file = paste0(file_path,"/2026-01-25_melanoma_filter.rds"))

meta_atlas8$major_cluster1 = case_when(
  meta_atlas8$RNA_snn_res.0.2 %in% c(8,9,10) ~ "caf-endo",
  meta_atlas8$RNA_snn_res.0.2 %in% c(1,4,3,7,5,6,12) ~ "melanoma cell",
  meta_atlas8$RNA_snn_res.0.2 %in% c(2,11) ~ "b-macro",
  meta_atlas8$RNA_snn_res.0.2 %in% c(0) ~ "t-nk"
)

meta_atlas8@meta.data$mega_merge = NA
meta_atlas8@meta.data$mega_merge <- apply(meta_atlas8@meta.data %>% select(major_cluster2, major_cluster1), 1, function(x) {
  if (max(table(x)) > 1) {
    # Get the most common annotation
    most_common <- names(which.max(table(x)))
    return(most_common)
  } else {
    return("not defined")
  }
})

meta_melanoma = FindClusters(meta_melanoma, resolution = 0.2)
meta_melanoma <- AddModuleScore(
  object = meta_melanoma,
  features = list(c("CD8A", "CD8B", "GZMB", "GNLY", "NKG7", "KLRB1")),
  name = "TNK_SigScore"
)

meta_atlas8 <- AddModuleScore(
  object = meta_atlas8,
  features = list(c(
    "LOC646214", "LOC643406", "LOC100131257", "LOC100190986", "LOC286437",
    "LOC100506190", "LOC284454", "LOC284023", "LOC100289019", "LOC100130231",
    "LOC100506385", "LOC100128288", "LOC100128531", "LOC100130581",
    "LOC100505876", "LOC729987", "LOC100506472"
  )),
  name = "Trash_Score"
)

meta_melanoma <- AddModuleScore(
  object = meta_melanoma,
  features = list(c("CD3D", "CD3E"
  )),
  name = "T_Score"
)

VlnPlot(meta_atlas8, features = "Trash_Score1") + 
  scale_y_continuous(breaks = seq(-0.5, 1, by = 0.05))

marker_list <- list(
  "Melanoma"    = c("PMEL", "MLANA", "MITF"),
  "T-CD8"       = c("CD8A", "CD8B", "GZMB"),
  "T-CD4"       = c("CD3D", "IL7R", "FOXP3","CCR7"),
  "NK"          = c("GNLY", "NKG7", "KLRB1"),
  "B-cell"      = c("MS4A1", "CD79A", "CD19"),
  "Macrophage"  = c("CD68", "CD163", "AIF1"),
  "Endothelial" = c("VWF", "PECAM1", "PLVAP"),
  "Mural−cell"  = c("RGS5", "ACTA2", "MCAM"),
  "CAF"         = c("COL1A1", "DCN", "FAP")
)

# Flatten the list into a simple vector of genes for Seurat
genes_to_plot <- unique(unlist(marker_list))

# 3. Generate the organized DotPlot
DotPlot(meta_melanoma, features = genes_to_plot, group.by = "RNA_snn_res.0.2") + 
  RotatedAxis() +                       # Rotate x-axis labels
  coord_flip() +                        # OPTIONAL: Flip to put Cell Types on X-axis, Genes on Y
  #scale_color_distiller(palette = "RdBu") +   # Nice Red-Blue color scheme
  ggtitle("Cell Type Markers") +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))
meta_melanoma <- AddModuleScore(
  object = meta_melanoma,
  features = list(c("CD3D", "CD3E"
  )),
  name = "T_Score"
)
p_meta <- DimPlot(
  meta_atlas9, 
  reduction = "umap", 
  group.by = c("RNA_snn_res.0.2", "Study", "mega_cluster"), 
  
  # This puts 10X_Public and SmartSeq2 ON TOP.
  # GSE215121_10X will be at the bottom by default.
  order = c("10X_Public", "GSE115978_SmartSeq2"), 
  
  raster = FALSE, 
  label = TRUE
) + NoLegend()

meta_fibend_markers = find_target_markers(meta_atlas8, c(12))
meta_fibend_markers = find_markers(meta_fibend5)


write.csv(meta_fibend_markers, paste0(file_path,"/melanoma_clusters.csv"), row.names = FALSE)


### Refine-t-nk cluster ----

meta_tnk = subset(meta_atlas9, RNA_snn_res.0.2 %in% c(6,7))
meta_tnk = recluster(meta_tnk)
meta_tnk = FindClusters(meta_tnk, resolution = 0.2)

meta_tnk2 = subset(meta_tnk, B_Score1 < 0.5 & Macrophage_Score1 < 0.5
                   & Melanoma_Score1 < 0.5 & Trash_Score1 < 0.5)
meta_tnk2 = recluster(meta_tnk2)
meta_tnk2 = FindClusters(meta_tnk2, resolution = 1.6)

meta_tnk3 = subset(meta_tnk2, RNA_snn_res.1.6 != 4)
meta_tnk3 = recluster(meta_tnk3)
meta_tnk3 = FindClusters(meta_tnk3, resolution = 2)

meta_tnk4 = subset(meta_tnk3, RNA_snn_res.2 %in% c(16,17,18,19), invert = T)
meta_tnk4 = recluster(meta_tnk4)
meta_tnk4 = FindClusters(meta_tnk4, resolution = 2)

meta_tnk5 = subset(meta_tnk4, RNA_snn_res.2 %in% c(16,17), invert = T)
meta_tnk5 = recluster(meta_tnk5)
meta_tnk5 = FindClusters(meta_tnk5, resolution = 2)

meta_tnk6 = subset(meta_tnk5, RNA_snn_res.2 %in% c(11,15), invert = T)
meta_tnk6 = recluster(meta_tnk6)
meta_tnk6 = FindClusters(meta_tnk6, resolution = 1)

meta_fibend_markers = find_target_markers(meta_tnk2, c(4))
meta_fibend_markers = find_markers(meta_tnk5)
write.csv(meta_fibend_markers %>% filter(p_val_adj <= 0.05), 
          paste0(file_path,"/melanoma_clusters.csv"), row.names = FALSE)
meta_tnk <- AddModuleScore(
  object = meta_tnk,
  features = list(c("CD68", "CD163", "AIF1")),
  name = "Macrophage_Score"
)
VlnPlot(meta_tnk, features = "Macrophage_Score1") + 
  scale_y_continuous(breaks = seq(-0.5, 1, by = 0.1))

meta_tnk <- AddModuleScore(
  object = meta_tnk,
  features = list(c("MS4A1", "CD79A", "CD19")),
  name = "B_Score"
)
VlnPlot(meta_tnk, features = "B_Score1") + 
  scale_y_continuous(breaks = seq(-0.5, 1, by = 0.1))

meta_tnk <- AddModuleScore(
  object = meta_tnk,
  features = list(c("MLANA","PMEL", "MLANA", "S100B")),
  name = "Melanoma_Score"
)
VlnPlot(meta_tnk, features = "Melanoma_Score1") + 
  scale_y_continuous(breaks = seq(-0.5, 1, by = 0.1))

meta_tnk <- AddModuleScore(
  object = meta_tnk,
  features = list(c(
    "LOC646214", "LOC643406", "LOC100131257", "LOC100190986", "LOC286437",
    "LOC100506190", "LOC284454", "LOC284023", "LOC100289019", "LOC100130231",
    "LOC100506385", "LOC100128288", "LOC100128531", "LOC100130581",
    "LOC100505876", "LOC729987", "LOC100506472"
  )),
  name = "Trash_Score"
)
VlnPlot(meta_tnk, features = "Trash_Score1") + 
  scale_y_continuous(breaks = seq(-0.5, 1, by = 0.1))

DotPlot(meta_tnk6, features = genes_to_plot, group.by = "RNA_snn_res.2") + 
  RotatedAxis() +                       # Rotate x-axis labels
  coord_flip() +                        # OPTIONAL: Flip to put Cell Types on X-axis, Genes on Y
  #scale_color_distiller(palette = "RdBu") +   # Nice Red-Blue color scheme
  ggtitle("Cell Type Markers") +
  theme(axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))
p_meta_tnk <- DimPlot(
  meta_tnk6, 
  reduction = "umap", 
  group.by = c("Study","RNA_snn_res.2", "mega_cluster"), 
  
  # This puts 10X_Public and SmartSeq2 ON TOP.
  # GSE215121_10X will be at the bottom by default.
  order = c("10X_Public", "GSE115978_SmartSeq2"), 
  
  raster = FALSE, 
  label = TRUE
) + NoLegend()

tnk_mela_diff = setdiff(colnames(meta_tnk),colnames(meta_tnk6))
meta_atlas10 = subset(meta_atlas9, cells = tnk_mela_diff, invert = T)
meta_atlas10 = recluster(meta_atlas10)
meta_atlas10 = FindClusters(meta_atlas10, resolution = 0.2)

NK_cells <- colnames(meta_tnk6)[meta_tnk6$RNA_snn_res.2 %in% c(1, 12, 4, 8)]
T_CD8s <- colnames(meta_tnk6)[meta_tnk6$RNA_snn_res.2 %in% c(2,11,13)]
T_CD4s <- colnames(meta_tnk6)[meta_tnk6$RNA_snn_res.2 %in% c(6,3,0,7,9,5,10)]
# 2. Assign "NK" to those cells in the main object
# Ensure the column exists and is a character vector (to avoid factor level errors)
meta_atlas10$mega_cluster[colnames(meta_atlas10) %in% NK_cells] <- "NK"
meta_atlas10$mega_cluster[colnames(meta_atlas10) %in% T_CD8s] <- "T-CD8"
meta_atlas10$mega_cluster[colnames(meta_atlas10) %in% T_CD4s] <- "T-CD4"


p_meta <- DimPlot(
  meta_atlas10, 
  reduction = "umap", 
  group.by = c("Study","RNA_snn_res.0.2", "mega_cluster"), 
  
  # This puts 10X_Public and SmartSeq2 ON TOP.
  # GSE215121_10X will be at the bottom by default.
  order = c("10X_Public", "GSE115978_SmartSeq2"), 
  
  raster = FALSE, 
  label = TRUE
) + NoLegend()

saveRDS(meta_atlas10, file = paste0(file_path,"/2026-01-25_melanoma_filter.rds"))

meta_atlas10$major_cluster1 = case_when(
  meta_atlas10$RNA_snn_res.0.2 %in% c(10,11,12) ~ "caf-endo",
  meta_atlas10$RNA_snn_res.0.2 %in% c(6,9,4,5,2,0,3) ~ "melanoma cell",
  meta_atlas10$RNA_snn_res.0.2 %in% c(1,13) ~ "b-macro",
  meta_atlas10$RNA_snn_res.0.2 %in% c(7,8) ~ "t-nk",
  meta_atlas10$RNA_snn_res.0.2 %in% c(14) ~ "not defined"
)

meta_atlas10$major_cluster2 <- case_when(
  meta_atlas10$mega_cluster %in% c("B_cell", "Macrophage") ~ "b-macro",
  meta_atlas10$mega_cluster %in% c("T-CD4", "T-CD8", "NK") ~ "t-nk",
  meta_atlas10$mega_cluster %in% c("Melanoma") ~ "melanoma cell",
  meta_atlas10$mega_cluster %in% c("CAF", "Mural_cell", "Endothelial") ~ "caf-endo"
)

meta_atlas10@meta.data$mega_merge = NA
meta_atlas10@meta.data$mega_merge <- apply(meta_atlas10@meta.data %>% select(major_cluster2, major_cluster1), 1, function(x) {
  if (max(table(x)) > 1) {
    # Get the most common annotation
    most_common <- names(which.max(table(x)))
    return(most_common)
  } else {
    return("not defined")
  }
})


meta_atlas11 = subset(meta_atlas10, mega_merge != "not defined")
meta_atlas11 = recluster(meta_atlas11)
meta_atlas11 = FindClusters(meta_atlas11, resolution = 0.2)

p_meta <- DimPlot(
  meta_atlas11, 
  reduction = "umap", 
  group.by = c("Study","RNA_snn_res.0.2", "mega_cluster"), 
  
  # This puts 10X_Public and SmartSeq2 ON TOP.
  # GSE215121_10X will be at the bottom by default.
  order = c("10X_Public", "GSE115978_SmartSeq2"), 
  
  raster = FALSE, 
  label = TRUE
) + NoLegend()

meta_atlas11$major_cluster1 = case_when(
  meta_atlas11$RNA_snn_res.0.2 %in% c(9) ~ "Endothelial",
  meta_atlas11$RNA_snn_res.0.2 %in% c(10) ~ "CAF",
  meta_atlas11$RNA_snn_res.0.2 %in% c(11) ~ "Mural_cell",
  meta_atlas11$RNA_snn_res.0.2 %in% c(4,6,5,2,0,3) ~ "Melanoma",
  meta_atlas11$RNA_snn_res.0.2 %in% c(1) ~ "B_cell",
  meta_atlas11$RNA_snn_res.0.2 %in% c(12) ~ "Macrophage",
  meta_atlas11$RNA_snn_res.0.2 %in% c(7,8) ~ "t-nk"
)

meta_atlas11$major_cluster2 <- case_when(
  meta_atlas11$mega_cluster %in% c("B_cell") ~ "B_cell",
  meta_atlas11$mega_cluster %in% c("Macrophage") ~ "Macrophage",
  meta_atlas11$mega_cluster %in% c("T-CD4", "T-CD8", "NK") ~ "t-nk",
  meta_atlas11$mega_cluster %in% c("Melanoma") ~ "Melanoma",
  meta_atlas11$mega_cluster %in% c("CAF") ~ "CAF",
  meta_atlas11$mega_cluster %in% c("Mural_cell") ~ "Mural_cell",
  meta_atlas11$mega_cluster %in% c("Endothelial") ~ "Endothelial"
)

meta_atlas11@meta.data$mega_merge = NA
meta_atlas11@meta.data$mega_merge <- apply(meta_atlas11@meta.data %>% select(major_cluster2, major_cluster1), 1, function(x) {
  if (max(table(x)) > 1) {
    # Get the most common annotation
    most_common <- names(which.max(table(x)))
    return(most_common)
  } else {
    return("not defined")
  }
})


meta_atlas12 = subset(meta_atlas11, mega_merge != "not defined")
meta_atlas12 = recluster(meta_atlas12)
meta_atlas12 = FindClusters(meta_atlas12, resolution = 0.2)
p_meta <- DimPlot(
  meta_atlas12, 
  reduction = "umap", 
  group.by = c("Study","RNA_snn_res.0.2", "mega_cluster"), 
  
  # This puts 10X_Public and SmartSeq2 ON TOP.
  # GSE215121_10X will be at the bottom by default.
  order = c("10X_Public", "GSE115978_SmartSeq2", "T-CD8"), 
  
  raster = FALSE, 
  label = TRUE
) + NoLegend()
saveRDS(meta_atlas12, file = paste0(file_path,"/2026-01-25_melanoma_filter.rds"))
### Make Pseudobulk and expression matrix ----
setwd(file_path)
meta_atlas12 = readRDS("../data/2026-01-25_melanoma_filter.rds")

melanoma_final_counts <- GetAssayData(object = meta_atlas12, assay = "RNA", layer = "counts")
melanoma_final_indent = data.frame(annotation = meta_atlas12@meta.data$mega_cluster, row.names = Cells(meta_atlas12))

# Normalize the cell numbers to percentages
melanoma_final_indent_sum <- data.frame(table(melanoma_final_indent$annotation))
colnames(melanoma_final_indent_sum) = c("annotation", "Freq")
melanoma_final_indent_sum$frac = melanoma_final_indent_sum$Freq / sum(melanoma_final_indent_sum$Freq)
sum(melanoma_final_indent_sum$Freq)

# Define the list of cell types for pseudobulk generation
cell_types <- unique(meta_atlas12@meta.data$mega_cluster)
cell_types = melanoma_final_indent_sum$annotation
# Set up the fraction list for sampling
fractionlist <- seq(0, 0.4, by = 0.02)

# Loop over each cell type for pseudobulk
for (cluster_cell_name in cell_types) {
  print(cluster_cell_name)
  # Locate the cells corresponding to the current cell type
  cluster_cell_coordinate <- which(melanoma_final_indent$annotation == cluster_cell_name)
  cluster_cells_barcodes <- rownames(melanoma_final_indent %>% filter(annotation == cluster_cell_name))
  
  # Remove these cells from the expression matrix
  melanoma_final_counts_remain <- melanoma_final_counts[, -which(colnames(melanoma_final_counts) %in% cluster_cells_barcodes)]
  num_remain_cells = ncol(melanoma_final_counts_remain)
  melanoma_final_counts_remain = rowSums(melanoma_final_counts_remain)
  
  # Initialize the pseudobulk expression matrix for the current cell type
  bulkexp_allfractions <- data.frame(row.names = row.names(melanoma_final_counts))
  bulkexp_allfractions$sample_control <- rowSums(melanoma_final_counts)
  
  # Perform bootstrap sampling for each fraction
  for (fraction in fractionlist) {
    
    num_of_cluster_cells = round(fraction * num_remain_cells / (1 - fraction))
    
    # Bootstrap the cluster cells for 10 replicates
    bootstrap_cells_count = data.frame(row.names = row.names(melanoma_final_counts))
    for (seed in seq(1, 10, by = 1)) {
      set.seed(seed)
      bootstrap_cells <- sample(cluster_cell_coordinate, num_of_cluster_cells, replace = TRUE)
      bootstrap_cells_count[, seed] <- rowSums(melanoma_final_counts[, bootstrap_cells])
    }
    melanoma_final_counts_cluster = rowSums(bootstrap_cells_count) / 10
    
    # Add pseudobulk sample to the expression matrix
    bulkexp_inloop <- rowSums(cbind(melanoma_final_counts_remain, melanoma_final_counts_cluster))
    cluster_colname = paste0(gsub(" ", "_", cluster_cell_name), gsub("0\\.", "", sprintf("%.3f", fraction)))
    bulkexp_allfractions[, cluster_colname] <- bulkexp_inloop
  }
  
  # Convert counts to CPM
  for (col in colnames(bulkexp_allfractions)) {
    newsum <- sum(bulkexp_allfractions[, col])
    bulkexp_allfractions[, col] <- (1000000 / newsum) * bulkexp_allfractions[, col]
  }
  
  # Add gene names and save to a file
  bulkexp_cpm_gene <- cbind(gene = row.names(bulkexp_allfractions), bulkexp_allfractions)
  output_filename <- paste0("../data/2026-01-28_pseudobulk_melanoma_", cluster_cell_name, "_cpm.txt")
  write.table(bulkexp_cpm_gene, output_filename, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
}


# Extract gene expression data for 4000 features
meta_atlas12 <- FindVariableFeatures(meta_atlas12, nfeatures = 4000)

coding_gene = read.table("../data/protein_coding_genes.txt", sep = "\t")$V1
var_genes <- VariableFeatures(meta_atlas12, nfeatures = 4000)

melanoma_obj_defined_feature <- data.frame(melanoma_final_counts[intersect(var_genes, coding_gene),], check.names = FALSE)
sum(row.names(melanoma_final_indent) == colnames(melanoma_obj_defined_feature))


# Change the col names
colnames(melanoma_obj_defined_feature) = melanoma_final_indent$annotation
#colnames(melanoma_obj_filtered_feature) <- gsub(" ", "_", colnames(melanoma_obj_filtered_feature))
unique(colnames(melanoma_obj_defined_feature))

# Write to a txt file
library(readr)
write_tsv(cbind(gene = rownames(melanoma_obj_defined_feature), melanoma_obj_defined_feature),"../data/2026-01-28_melanoma_obj_feature_4000.txt")


### Analyze pseudo bulk result ----
library(DescTools)
library(tidyr)
all_files <- c(list.files("../data/2026-01-28_melanoma_pseudobulk", full.names = TRUE, recursive = TRUE))
pseudo_result_files <- all_files[grepl("CIBERSORTx", all_files)]

pseudo_summary <- data.frame(
  Cell = character(),
  Slope = numeric(),
  P_Value_Slope = numeric(),
  R_squared = numeric(),
  Spearman_Rho = numeric(), # Added for Spearman correlation
  RMSE = numeric(),         # Added for Root Mean Square Error
  CCC = numeric(),          # Added for Lin's Concordance Correlation Coefficient
  stringsAsFactors = FALSE
)
for (file in pseudo_result_files) {
  cell_name <- sub(".*_800_(.*?)_pseudo.*", "\\1", file)
  pseudo_result = read.table(file, sep = "\t", header = T, check.names = F)
  pseudo_result = pseudo_result[-1, ]
  pseudo_result$truth <- substr(pseudo_result$Mixture, nchar(pseudo_result$Mixture)-2, nchar(pseudo_result$Mixture))
  pseudo_result$truth <- as.numeric(pseudo_result$truth)/1000
  plot = ggplot(pseudo_result, aes(y = .data[[cell_name]], x = truth)) + geom_point() + 
    labs(
      title = cell_name,  # Add a meaningful title
      x = "Truth",                                  # Label for x-axis
      y = "Cibersort"                              # Label for y-axis (using the cell name dynamically)
    )
  ggsave(paste0("../data/2026-01-28_melanoma_pseudobulk/", cell_name, ".pdf"), plot)
  # Perform linear regression using the actual column reference
  model <- lm(pseudo_result[[cell_name]] ~ pseudo_result$truth)
  slope <- coef(model)[2]
  p_value_slope <- summary(model)$coefficients[2, 4]
  r_squared <- summary(model)$r.squared
  
  # 2. Spearman Correlation
  spearman_test <- cor.test(pseudo_result[[cell_name]], pseudo_result$truth, method = "spearman", use = "complete.obs")
  spearman_rho <- spearman_test$estimate
  
  # 3. Root Mean Square Error (RMSE)
  rmse <- sqrt(mean((pseudo_result[[cell_name]] - pseudo_result$truth)^2, na.rm = TRUE))
  
  # 4. Lin's Concordance Correlation Coefficient (CCC)
  #    Requires the 'DescTools' package
  ccc_result <- CCC(pseudo_result[[cell_name]], pseudo_result$truth)
  ccc_value <- ccc_result$rho.c$est # Extract the estimate
  
  # --- Append all results to the summary data frame ---
  pseudo_summary <- rbind(pseudo_summary, data.frame(
    Cell = cell_name,
    Slope = slope,
    P_Value_Slope = p_value_slope, # Renamed for clarity
    R_squared = r_squared,
    Spearman_Rho = spearman_rho,
    RMSE = rmse,
    CCC = ccc_value
  ))
}

# Create a new column for -log10(P_value)
pseudo_summary$log_p_value <- -log10(pseudo_summary$P_Value)

# Generate the scatter plot
p_pseudo = ggplot(pseudo_summary, aes(x = Cell, y = Slope, color = log_p_value)) +
  geom_point(size = 4) +  # Increase the point size for better visibility
  scale_color_gradient(low = "blue", high = "red") +  # Color gradient for p-values
  labs(
    title = "Slope vs Cell Type with -log10(P_value)",
    x = "Cell Type",
    y = "Slope",
    color = "-log10(P_value)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rota

p_pseudo = ggplot(pseudo_summary%>%filter(Filter == "Filter"), aes(x = Cell, y = Slope, color = log_p_value)) +
  geom_point(size = 4) +  # Increase the point size for better visibility
  scale_color_gradient(low = "blue", high = "red") +  # Color gradient for p-values
  labs(
    title = "Slope vs Cell Type with -log10(P_value)",
    x = "Cell Type",
    y = "Slope",
    color = "-log10(P_value)"
  ) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rota

pseudo_summary_long <- pseudo_summary %>%
  select(Cell, Slope, R_squared, Spearman_Rho, CCC, RMSE) %>%
  pivot_longer(
    cols = -Cell, 
    names_to = "Metric", 
    values_to = "Value"
  ) %>%
  # Make Metric a factor to control the order of facets in the plot
  mutate(Metric = factor(Metric, levels = c("Slope", "R_squared", "Spearman_Rho", "CCC", "RMSE")))

# Create the faceted bar plot
performance_plot <- ggplot(pseudo_summary_long, aes(x = reorder(Cell, Value), y = Value, fill = Metric)) +
  geom_col() + # Use geom_col for bar chart where height is the value
  facet_wrap(~ Metric, scales = "free_y") + # Create separate panels for each metric with their own y-axis scale
  coord_flip() + # Flip coordinates to make cell names easier to read
  labs(
    title = "Deconvolution Performance Metrics for melanoma pseudobulk",
    x = "Cell Type",
    y = "Value"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position = "none", # The fill color is redundant with the facet title
    strip.background = element_rect(fill = "grey85"),
    strip.text = element_text(face = "bold")
  )

print(performance_plot)


### Make fraction plot ----
metadata = read.table("../data/2025/2025-01-27_GeoMxmetadata.txt", sep = "\t", header = T)
fraction = read.table("../data/2026-01-28_cibersort_fraction.txt", sep = "\t", header = T)


cell_name = colnames(fraction)[2:10]
frac_meta = merge(metadata, fraction, by.x = "ID", by.y = "Mixture")

# Reorder levels first, if needed
frac_meta$Pathophysiology <- factor(
  frac_meta$Pathophysiology,
  levels = c("Norm", "Nevi", "Melanoma")
)


# Create a list of pairwise comparisons
my_comparisons <- list(
  c("Norm", "Nevi"),
  c("Nevi", "Melanoma"),
  c("Norm", "Melanoma")
)

for (cell in cell_name) {
  
  p <- ggplot(frac_meta, aes(
    x = Pathophysiology,
    y = .data[[cell]],
    fill = Segment_Tags
  )) +
    geom_boxplot(
      alpha = 0.7,
      outlier.shape = NA,
      position = position_dodge(width = 0.8)
    ) +
    geom_jitter(
      position = position_jitterdodge(
        dodge.width = 0.8,
        jitter.width = 0.2
      ),
      alpha = 0.8,
      size = 1.5
    ) +
    theme_bw(base_size = 12) +
    labs(
      title = paste0("Comparison of ", cell, " Fractions"),
      x = "Pathophysiology",
      y = paste0(cell, " Fraction")
    ) +
    scale_fill_discrete(name = "Pathophysiology") +
    
    # 1) Overall p-value (Kruskal-Wallis) by default since >2 groups:
    stat_compare_means(
      label = "p.format",
      label.y = max(frac_meta[[cell]], na.rm = TRUE) * 1.1,  # Position p-value above the top
      size = 4
    ) #+
  
  
    
    # 2) Pairwise comparisons using Wilcoxon test:
    #stat_compare_means(
    #  comparisons = my_comparisons,
    #  method = "wilcox.test",
    #  label = "p.signif",        # Show significance (e.g., "***", "**", etc.)
    #  # label = "p.format",      # Alternatively show numeric p-values
    #  tip.length = 0.01,
    #  vjust = 0.5
    #)
  
  # Save the plot
  ggsave(
    filename = paste0("../analysis/2026-01-28_patho_", cell, "_Fractions.pdf"),
    plot = p,
    width = 6,   # Adjust the width as needed
    height = 4,  # Adjust the height as needed
    units = "in",  # units can be "in", "cm", or "mm"
    dpi = 300      # resolution
  )
}

frac_meta$Path_sub_type <- factor(
  frac_meta$Path_sub_type,
  levels = c("Norm", "Comp", "ID","Junc", "AMEL", "LMM", "MET", "NOD", "SITU", "SSM")
)

for (cell in cell_name) {
  
  p <- ggplot(frac_meta, aes(
    x = Path_sub_type,
    y = .data[[cell]],
    fill = Segment_Tags
  )) +
    geom_boxplot(
      alpha = 0.7,
      outlier.shape = NA,
      position = position_dodge(width = 0.8)
    ) +
    geom_jitter(
      position = position_jitterdodge(
        dodge.width = 0.8,
        jitter.width = 0.2
      ),
      alpha = 0.8,
      size = 1.5
    ) +
    theme_bw(base_size = 12) +
    labs(
      title = paste0("Comparison of ", cell, " Fractions"),
      x = "Path_sub_type",
      y = paste0(cell, " Fraction")
    ) +
    scale_fill_discrete(name = "Path_sub_type") +
    
    # 1) Overall p-value (Kruskal-Wallis) by default since >2 groups:
    stat_compare_means(
      label = "p.format",
      label.y = max(frac_meta[[cell]], na.rm = TRUE) * 1.1,  # Position p-value above the top
      size = 4
    ) #+
  
  # 2) Pairwise comparisons using Wilcoxon test:
  #stat_compare_means(
  #  comparisons = my_comparisons,
  #  method = "wilcox.test",
  #  label = "p.signif",        # Show significance (e.g., "***", "**", etc.)
  #  # label = "p.format",      # Alternatively show numeric p-values
  #  tip.length = 0.01,
  #  vjust = 0.5
  #)
  
  # Save the plot
  ggsave(
    filename = paste0("../analysis/2026-01-28_patho_subtype_", cell, "_Fractions.pdf"),
    plot = p,
    width = 12,   # Adjust the width as needed
    height = 4,  # Adjust the height as needed
    units = "in",  # units can be "in", "cm", or "mm"
    dpi = 300      # resolution
  )
}
######## Make melanoma only data ----
metadata = read.table("../data/2025/2025-01-27_GeoMxmetadata.txt", sep = "\t", header = T)

melig_id = metadata %>% filter(Path_type != "Norm")
exp = read.table("../data/2025/2025-01-24_GeoMx_expression_data.txt", sep = "\t", header = T)
exp_noNorm = exp[,c("Gene", melig_id$ID)]
write.table(exp_noNorm, "2025/2025-01-24_GeoMx_expression_data_noNorm.txt", sep = "\t", quote = F, row.names = F)


######## Analyze melanoma hires data ----
metadata = read.table("../data/2025/2025-01-27_GeoMxmetadata.txt", sep = "\t", header = T) %>% 
  filter(Path_type != "Norm", Segment_Tags == "S100") %>% select(ID, Path_type, Path_sub_type)
melanoma_rec = read.table("../data/2026-01-28_GeoMx_recovery/Melanoma/heatmap_data.txt", sep = "\t", header = T) %>%
  filter(apply(., 1, var) != 0)

melanoma_hires = read.table("../data/2026-01-28_melanoma_GeoMx_hires/CIBERSORTxHiRes_NA_Melanoma_Window36.txt", 
                            sep = "\t", header = T, row.names = 1) %>%
  filter(apply(., 1, var) != 0)

# Load the correct packages for GeoMx / Bulk data
library(limma)
library(dplyr)
library(pheatmap)

# ---------------------------------------------------------
# 1. Align and format the metadata
# ---------------------------------------------------------
rownames(metadata) <- metadata$ID

# Force Path_type to be a factor with your exact desired order
metadata$Path_type <- factor(metadata$Path_type, levels = c("Nevi", "MALIGN", "MET"))

# Sort the metadata dataframe by Path_type first, then Path_sub_type
metadata_ordered <- metadata[order(metadata$Path_type, metadata$Path_sub_type), ]

# Reorder the expression matrix columns to match this perfectly sorted metadata
# (Assuming your matrix from earlier is called 'melanoma_hires' or 'exp_noNorm')
expr_mat_ordered <- melanoma_hires[, rownames(metadata_ordered)]

# GeoMx data usually needs to be log2-transformed before running linear models
# (If your data is already log2 transformed, skip the '+ 1' math here)
log_expr <- log2(expr_mat_ordered + 1)

# ---------------------------------------------------------
# 2. GeoMx-Specific Differential Expression (One-vs-All loop)
# ---------------------------------------------------------
subtypes <- unique(metadata_ordered$Path_sub_type)
all_markers <- data.frame()

for(st in subtypes) {
  # Create a binary comparison: Current Sub-type vs Everything Else
  group <- ifelse(metadata_ordered$Path_sub_type == st, "Target", "Other")
  group <- factor(group, levels = c("Other", "Target"))
  
  # Build the linear model
  design <- model.matrix(~ group)
  fit <- lmFit(log_expr, design)
  fit <- eBayes(fit) # Applies empirical Bayes smoothing (standard for limma)
  
  # Extract the results for Target vs Other
  res <- topTable(fit, coef = "groupTarget", number = Inf)
  res$gene <- rownames(res)
  res$Path_sub_type <- st
  
  # Keep only genes that are significantly UPREGULATED in this sub-type
  res_pos <- res %>% filter(logFC > 0.25, adj.P.Val < 0.05)
  
  all_markers <- rbind(all_markers, res_pos)
}

# ---------------------------------------------------------
# 3. Extract the Top 100 DEGs
# ---------------------------------------------------------
# Pull the top ~11 genes per sub-type to get exactly ~100 distinct markers total
top_genes <- all_markers %>%
  group_by(Path_sub_type) %>%
  slice_max(n = 11, order_by = logFC) %>%
  pull(gene) %>%
  unique()

# ---------------------------------------------------------
# 4. Prepare and Plot the Heatmap
# ---------------------------------------------------------
# Extract only the top genes from the log2 matrix
heatmap_matrix <- log_expr[top_genes, ]

# Create the annotation dataframe for the top of the heatmap
anno_col_df <- data.frame(
  Path_type = metadata_ordered$Path_type,
  Path_sub_type = metadata_ordered$Path_sub_type,
  row.names = rownames(metadata_ordered)
)

sub_type_counts <- rle(as.character(metadata_ordered$Path_sub_type))$lengths

# cumsum() turns those counts into exact column indices for the breaks
col_gap_indices <- cumsum(sub_type_counts)

# Remove the very last index (pheatmap will throw an error if you tell it to put a gap after the final column)
col_gap_indices <- col_gap_indices[-length(col_gap_indices)]
# ------------------------------------------------

# Plot the carefully ordered heatmap
pheatmap(
  mat = as.matrix(heatmap_matrix),
  scale = "row",                  
  cluster_cols = FALSE,           # CRITICAL: gaps_col ONLY works when cluster_cols = FALSE
  cluster_rows = TRUE,            
  annotation_col = anno_col_df,
  gaps_col = col_gap_indices,     # <--- Adds the visual gaps between sub-types
  show_colnames = FALSE,          
  fontsize_row = 6,               
  main = "Top GeoMx DEGs Across Melanoma Progression"
)


## PCA Plot

# Load ggplot2 for visualization
library(ggplot2)

# ---------------------------------------------------------
# 1. Prepare the Data for PCA
# ---------------------------------------------------------
# Subset the log-transformed expression matrix to ONLY include our DEG list
# (Using the 'top_genes' vector and 'log_expr' matrix from the previous limma step)


deg_expr <- log_expr[unique(all_markers$gene), ]

# prcomp() expects samples as rows and genes as columns, so we must transpose (t) the matrix
pca_input <- t(deg_expr)

# Run the Principal Component Analysis
# scale. = TRUE ensures genes with massive expression don't overpower lower-expressed genes
pca_res <- prcomp(pca_input, scale. = TRUE)

# ---------------------------------------------------------
# 2. Extract Data and Merge with Metadata
# ---------------------------------------------------------
# Extract the principal component coordinates for each sample
pca_df <- as.data.frame(pca_res$x)

# Safely add the metadata back in by matching rownames
pca_df$ID <- rownames(pca_df)
pca_df <- merge(pca_df, metadata_ordered, by = "ID")

# Calculate the percentage of variance explained by PC1 and PC2 (for the axis labels)
var_explained <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)

library(dplyr)
library(ggplot2)

# ---------------------------------------------------------
# 1. Calculate the center (centroid) of each Sub-type
# ---------------------------------------------------------
centroids <- pca_df %>%
  group_by(Path_sub_type) %>%
  summarise(
    PC1 = mean(PC1),
    PC2 = mean(PC2),
    # Keep the parent Path_type so the centroid shapes match the background dots
    Path_type = unique(Path_type)[1] 
  )

# ---------------------------------------------------------
# 2. Define the Biological Progression Trajectory
# ---------------------------------------------------------
# CRITICAL: R needs to know the order in which to draw the connecting line.
# Please adjust this list to match your exact biological hypothesis 
# (e.g., from early Nevi -> in situ -> invasive -> Metastasis)
progression_order <- c("Comp", "ID", "Junc", "SITU", "LMM", "SSM", "NOD", "AMEL", "MET")

# Apply this strict order to the centroids so the line flows correctly
centroids <- centroids %>%
  mutate(Path_sub_type = factor(Path_sub_type, levels = progression_order)) %>%
  arrange(Path_sub_type)

# ---------------------------------------------------------
# 3. Generate the Trajectory PCA Plot
# ---------------------------------------------------------
ggplot(pca_df, aes(x = PC1, y = PC2)) +
  # 1. Plot the individual patient ROIs (made slightly transparent to reduce clutter)
  geom_point(aes(color = Path_sub_type, shape = Path_type), size = 2.5, alpha = 0.3) +
  
  # 2. Plot the calculated centroids as larger, bold points
  geom_point(data = centroids, aes(color = Path_sub_type, shape = Path_type), 
             size = 6, stroke = 1.5) +
  
  # 3. Draw the directed progression path connecting the centroids!
  geom_path(data = centroids, aes(group = 1), color = "black", linewidth = 1,
            arrow = arrow(type = "closed", length = unit(0.15, "inches"))) +
  
  # Clean, publication-ready theme
  theme_classic() +
  labs(
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)"),
    title = "Melanoma Progression Trajectory in PCA Space",
    subtitle = "Arrows indicate the mean transcriptomic shift between biological stages",
    color = "Sub-type",
    shape = "Tumor Type"
  ) +
  
  # Text formatting
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30"),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )


# ---------------------------------------------------------
# 3. Generate the Plot without trajectory
# ---------------------------------------------------------
ggplot(pca_df, aes(x = PC1, y = PC2, color = Path_sub_type, shape = Path_type)) +
  # Add the sample dots
  geom_point(size = 3.5, alpha = 0.85) +
  
  # Optional: Draw a subtle outline around the shapes to make them pop
  # (Requires using shapes 21-25 if you want custom fill/borders, but default works great too)
  
  # Clean, publication-ready theme
  theme_classic() +
  
  # Custom axis labels showing exact variance percentages
  labs(
    x = paste0("PC1 (", var_explained[1], "%)"),
    y = paste0("PC2 (", var_explained[2], "%)"),
    title = "PCA of GeoMx Samples",
    subtitle = "Clustering based on Top Differentially Expressed Genes",
    color = "Sub-type",
    shape = "Tumor Type"
  ) +
  
  # Refine text sizes and legend placement
  theme(
    plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "grey30"),
    axis.title = element_text(face = "bold", size = 12),
    axis.text = element_text(size = 10),
    legend.position = "right",
    legend.title = element_text(face = "bold")
  )


######## Functions needed ----
recluster <- function(seurat_obj) {
  
  # 1. Normalize (FIXED: use seurat_obj)
  seurat_obj <- NormalizeData(seurat_obj)
  
  # 2. Find Variable Features
  seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 3000)
  
  # 3. Scale Data
  seurat_obj <- ScaleData(seurat_obj)
  
  # 4. Run PCA
  seurat_obj <- RunPCA(seurat_obj, npcs = 20, verbose = FALSE)
  
  # 5. Run Harmony
  # Ensure you have the 'Study' column in metadata
  seurat_obj <- RunHarmony(seurat_obj, group.by.vars = "Study", plot_convergence = FALSE)
  
  # 6. Find Neighbors (Build Graph)
  seurat_obj <- FindNeighbors(
    seurat_obj,
    reduction = "harmony",
    dims      = 1:15  # increased slightly to capture more variance
  )
  
  # 7. Run UMAP
  seurat_obj <- RunUMAP(
    seurat_obj,
    reduction = "harmony",
    dims      = 1:15
  )
  
  return(seurat_obj)
}

find_markers <- function(seurat_obj) {
  
  all_tumor_markers <- FindAllMarkers(object = seurat_obj, 
                                      only.pos = TRUE, 
                                      min.pct = 0.25, 
                                      logfc.threshold = 0.25)
  
  top100_tumor <- all_tumor_markers %>%
    group_by(cluster) %>%
    slice_min(n = 100, order_by = p_val_adj)
  
  return(top100_tumor)
}

find_target_markers <- function(seurat_obj, target_clusters) {
  
  # Loop through them and calculate markers one by one
  marker_list <- lapply(target_clusters, function(id) {
    print(paste("Processing Cluster", id))
    
    # Check if cluster exists in the object's current identities to prevent errors
    if (id %in% unique(Idents(seurat_obj))) {
      
      # Run FindMarkers
      markers <- FindMarkers(seurat_obj, 
                             ident.1 = id, 
                             only.pos = TRUE, 
                             min.pct = 0.25, 
                             logfc.threshold = 0.25)
      
      # Add columns to match FindAllMarkers format
      markers$cluster <- id
      markers$gene <- rownames(markers)
      return(markers)
      
    } else {
      print(paste("Skipping Cluster", id, "- not found in object"))
      return(NULL)
    }
  })
  
  # Filter out NULL entries (clusters that were skipped)
  marker_list <- marker_list[!sapply(marker_list, is.null)]
  
  # Check if we have any results to combine
  if (length(marker_list) == 0) {
    message("No markers found for the specified clusters.")
    return(NULL)
  }
  
  # Combine into a single dataframe
  all_markers <- do.call(rbind, marker_list)
  
  # Select top 100 markers per cluster by adjusted p-value
  top100 <- all_markers %>%
    group_by(cluster) %>%
    slice_min(n = 100, order_by = p_val_adj)
  
  return(top100)
}
