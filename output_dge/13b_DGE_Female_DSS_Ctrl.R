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
  library(scDblFinder)
  library(BiocSingular)
  library(Libra)
})

emilia_seurat <- readRDS("output_NoDRD/06_Final_noDRD_slim.rds")

head(emilia_seurat)

table(emilia_seurat$Genotype)

emilia_seurat@meta.data$Diff[emilia_seurat@meta.data$Genotype == "Ctrl"] <- "B_Ctrl"
emilia_seurat@meta.data$Diff[emilia_seurat@meta.data$Genotype == "DSS"] <- "A_DSS"

subFemale <- subset(emilia_seurat, subset = Eif2s3y_Sex_Prediction == "Female")
table(subFemale$Diff)

tmp <- subFemale

tmp@meta.data$cell_type <- tmp@meta.data$Cell
tmp@meta.data$replicate <- tmp@meta.data$Genotype
tmp@meta.data$label <- tmp@meta.data$Diff

dir.create("DGE_NoDRD")

DE_DSS_vs_Ctrl_Female <- run_de(tmp, de_family = 'singlecell', de_method = 'wilcox',n_threads = 20)

saveRDS(DE_DSS_vs_Ctrl_Female, file = "DGE_NoDRD/DE_DSS_vs_Ctrl_Female.rds")

tmp <- DE_DSS_vs_Ctrl_Female %>%
  mutate(Abs = abs(avg_logFC)) %>%
  dplyr::filter(p_val_adj < 0.05 & Abs > 0.0) %>%
  mutate(Direction = case_when(avg_logFC > 0.0 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.0 & p_val_adj < 0.05 ~ "DownReg")) %>%
  arrange(desc(Abs))

openxlsx::write.xlsx(tmp, 
                     file = "DGE_NoDRD/Female_01_DGE_NoDRD_DSS_vs_Ctrl_Full.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="DGE_NoDRD",
                     overwrite=T)

DE_DSS_vs_Ctrl_Female_sign <- DE_DSS_vs_Ctrl_Female %>%
  mutate(Abs = abs(avg_logFC)) %>%
  dplyr::filter(p_val_adj < 0.05 & Abs > 0.3) %>%
  mutate(Direction = case_when(avg_logFC > 0.3 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.3 & p_val_adj < 0.05 ~ "DownReg")) %>%
  arrange(desc(Abs))

openxlsx::write.xlsx(DE_DSS_vs_Ctrl_Female_sign, 
                     file = "DGE_NoDRD/Female_01_DGE_NoDRD_DSS_vs_Ctrl_Filtered_pvalLT05_AbsGT03.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="DGE_NoDRD",
                     overwrite=T)

DE_DSS_vs_Ctrl_Female_sign <- DE_DSS_vs_Ctrl_Female %>%
  mutate(Abs = abs(avg_logFC)) %>%
  dplyr::filter(p_val_adj < 0.05 & Abs > 0.2) %>%
  mutate(Direction = case_when(avg_logFC > 0.2 & p_val_adj < 0.05 ~ "UpReg", avg_logFC < -0.2 & p_val_adj < 0.05 ~ "DownReg")) %>%
  arrange(desc(Abs))

openxlsx::write.xlsx(DE_DSS_vs_Ctrl_Female_sign, 
                     file = "DGE_NoDRD/Female_01_DGE_NoDRD_DSS_vs_Ctrl_Filtered_pvalLT05_AbsGT02.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="DGE_NoDRD",
                     overwrite=T)
