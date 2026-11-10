suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(emmeans)
  library(SingleCellExperiment)
  library(scater)
  library(BiocParallel)
  library(ggpubr)
  library(speckle)
  library(magrittr)
  library(broom)
  library(muscat)
  library(Seurat)
  library(clustree)
  library(leiden)
  library(data.table)
  library(cowplot)
  library(CellChat)
  library(BiocSingular)
  library(ComplexHeatmap)
  library(Augur)
})

# Cell type prioritization in high-dimensional single-cell data
setwd("/Users/suganyasubramanian/Emilia")

emilia_seurat <- readRDS("output_NoDRD/06_Final_noDRD_slim.rds")

head(emilia_seurat)

dir.create("output_Augur")

augur = calculate_auc(emilia_seurat, 
                      cell_type_col = "Cell", label_col = "Genotype",
                      n_threads = 16)

save(augur,file = "output_Augur/Augur.RData")

pdf("output_Augur/UMAP_Augur_red.pdf", width = 3, height = 3)
plot_umap(
  augur,
  emilia_seurat,
  mode = c("default", "rank"),
  reduction = "umap",
  #palette = "cividis",
  palette = c("white","red"),
  top_n = 0,
  cell_type_col = "Cell"
)
dev.off()

#Try scaling
aug_scaled <- scale(augur$AUC$auc)
aug_scaled

augur$AUC$auc <- aug_scaled

pdf("output_Augur/UMAP_Augur_scaled.pdf", width = 3, height = 3)
plot_umap(
  augur,
  emilia_seurat,
  mode = c("default", "rank"),
  reduction = "umap",
  palette = "cividis",
  top_n = 0,
  cell_type_col = "Cell"
)
dev.off()

#with lr
emilia_seurat <- readRDS("output_NoDRD/06_Final_noDRD_slim.rds")

augur1 = calculate_auc(emilia_seurat,
                       cell_type_col = "Cell", 
                       label_col = "Genotype",
                      classifier = "lr", #logistic reg rather than random forest (rf)
                      n_threads = 10)

save(augur1,file = "output_Augur/Augur_rdfor.RData")

pdf("output_Augur/UMAP_Augur1.pdf", width = 3, height = 3)
plot_umap(augur1, emilia_seurat, mode = "rank",
          palette = c("white","red"),cell_type_col = "Cell")
dev.off()

pdf("output_Augur/UMAP_Augur1_lollipop.pdf", width = 3, height = 3)
plot_lollipop(augur1)
dev.off()

p1 <- plot_umap(augur1, emilia_seurat, mode = "rank",
          palette = c("white","red"),cell_type_col = "Cell")
p2 <- DimPlot(object = emilia_seurat, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p3 <- plot_lollipop(augur1)

pdf("output_Augur/Augur_Final_Result.pdf", width = 12, height = 3)
plot_grid(p2, p1, p3, ncol =3)
dev.off()

