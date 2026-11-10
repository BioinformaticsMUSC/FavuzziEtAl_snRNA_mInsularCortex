suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(emmeans)
  library(SingleCellExperiment)
  library(scater)
  library(BiocParallel)
  library(ggpubr)
  #library(speckle)
  #library(magrittr)
  #library(broom)
  #library(muscat)
  library(Seurat)
  library(clustree)
  library(leiden)
  library(data.table)
  library(cowplot)
  library(scDblFinder)
  library(BiocSingular)
  #library(scds)
})

setwd("/Users/SuganyaSubramanian/Emilia")

load("output_Relabel/04_Final_slim_with_ClassifySex.RData") #emilia_seurat

source("/Users/SuganyaSubramanian/Utils.R")

table(emilia_seurat$Cell)
emilia_noDRD <- subset(x = emilia_seurat, subset = Cell == "D1/D2_STR", invert = TRUE)
emilia_noDRD <- subset(x = emilia_noDRD, subset = Drd1 == 0 & Drd2 == 0)
emilia_noDRD

emilia_noDRD$Ori_cluster <- emilia_noDRD$seurat_clusters

emilia_noDRD <- processing_seurat_sctransform(emilia_noDRD, 
                                                       vars_to_regress = c("nCount_RNA","pMito","pRibo"), 
                                                       npcs = 30, 
                                                       res = c(0.3,0.4,0.5,0.6,0.7,0.8))

DefaultAssay(emilia_noDRD) <- "RNA"
emilia_noDRD <- NormalizeData(object = emilia_noDRD, 
                                       normalization.method = "LogNormalize", 
                                       scale.factor = 10000)

dir.create("output_NoDRD")

table(emilia_noDRD$Cell)
Idents(emilia_noDRD) <- "Cell"

emilia_noDRD <- RenameIdents(object = emilia_noDRD, `L6_CT_Foxp2` = "Inh_Adarb2")

table(Idents(emilia_noDRD))

emilia_noDRD@meta.data$Cell <- emilia_noDRD@active.ident

save(emilia_noDRD, file = "output_NoDRD/05_Final_noDRD.RData")

pdf("output_NoDRD/10_UMAP_NoDRD_res03.pdf", width = 10, height = 6)
p1 <- DimPlot(object = emilia_noDRD, group.by = "SCT_snn_res.0.3",reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = emilia_noDRD, group.by = "Cell",reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2)
dev.off()

pdf("output_NoDRD/10_UMAP_NoDRD_res04.pdf", width = 10, height = 6)
p1 <- DimPlot(object = emilia_noDRD, group.by = "SCT_snn_res.0.4", reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = emilia_noDRD, group.by = "Cell",reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2)
dev.off()

pdf("output_NoDRD/10_UMAP_NoDRD_res05.pdf", width = 10, height = 6)
p1 <- DimPlot(object = emilia_noDRD, group.by = "SCT_snn_res.0.5",reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = emilia_noDRD, group.by = "Cell",reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2)
dev.off()

pdf("output_NoDRD/10_UMAP_NoDRD_res06.pdf", width = 10, height = 6)
p1 <- DimPlot(object = emilia_noDRD, group.by = "SCT_snn_res.0.6",reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = emilia_noDRD, group.by = "Cell",reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2)
dev.off()

pdf("output_NoDRD/10_UMAP_NoDRD_res07.pdf", width = 10, height = 6)
p1 <- DimPlot(object = emilia_noDRD, group.by = "SCT_snn_res.0.7",reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = emilia_noDRD, group.by = "Cell",reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2)
dev.off()

pdf("output_NoDRD/10_UMAP_NoDRD_res08.pdf", width = 10, height = 6)
p1 <- DimPlot(object = emilia_noDRD, group.by = "SCT_snn_res.0.8",reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = emilia_noDRD, group.by = "Cell",reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2)
dev.off()

pdf("output_NoDRD/01_UMAP_Labelled.pdf", width = 7, height = 6)
DimPlot(object = emilia_noDRD, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
dev.off()

pdf("output_NoDRD/01_UMAP_Labelled_CellClass.pdf", width = 7, height = 6)
DimPlot(object = emilia_noDRD, reduction = "umap", label = TRUE, pt.size = 0.5, group.by="Cell_Class") + theme(legend.position="none")
dev.off()

pdf("output_NoDRD/02_UMAP_Labelled_Groupby.pdf", width = 14, height = 6)
p1 <- DimPlot(object = emilia_noDRD, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = emilia_noDRD, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

pdf("output_NoDRD/03_UMAP_Labelled_withclusters.pdf", width = 14, height = 6)
Idents(emilia_noDRD) <- "seurat_clusters"
p1 <- DimPlot(object = emilia_noDRD, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
Idents(emilia_noDRD) <- "Cell"
p2 <- DimPlot(object = emilia_noDRD, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2)
dev.off()

pdf("output_NoDRD/04_UMAP_Labelled_Splitby_Genotype.pdf", width = 14, height = 6)
DimPlot(object = emilia_noDRD, reduction = "umap", label = TRUE,ncol = 2, pt.size = 0.5,split.by = "Genotype") + theme(legend.position="none")
dev.off()

order <- c("Astro_Gja1","Oligo_Mbp","OPC_Vcan","Micro_Dock8",
           "L2/3_IT_Olfm3","L2/3_IT_Lrrtm4","L4/5_IT_Pcsk5","L4/5_IT_Rorb","L4/5_IT_Il1rapl2",
           ,"L6_CT_Hs3st4","L6b_Ndst4","L6b_ENT_Tshz2",#"D1/D2_STR","L6_CT_Foxp2"
           "Inh_Pvalb","Inh_Adarb2","Inh_Sst",
           "Endo_Cldn5","VLMC_Slc6a13")

Idents(emilia_noDRD) <- factor(Idents(emilia_noDRD), levels = order) # reorder the factor based on cell type

markers <- c("Gja1","Mbp","Vcan","Dock8",
             "Olfm3","Lrrtm4","Pcsk5","Rorb","Il1rapl2",
             "Foxp2","Hs3st4","Ndst4","Tshz2",#"Drd1","Drd2",
             "Pvalb","Adarb2","Sst",
             "Cldn5","Slc6a13")

pdf("output_NoDRD/07_DotPlot_Celltype_Markers.pdf", width = 10, height = 6)
DotPlot(emilia_noDRD, features = markers) + rotate_x_text(angle = 45)
dev.off()

pdf("output_NoDRD/07_DotPlot_seurat_Markers.pdf", width = 10, height = 6)
DotPlot(emilia_noDRD, features = markers, group.by = "SCT_snn_res.0.7") + rotate_x_text(angle = 45)
dev.off()

### 2 Cell proportion ####
emilia_noDRD$Cell <- factor(emilia_noDRD$Cell, levels = order) # reorder the factor based on cell type

color <- c("#17154f","#b0799a","#e69b00","#355828","#4f2400","#b62684")
#"#e69b00","#355828","#5eb0e5","#ee7762","#4f2400","#b62684")

prop_cell <-  emilia_noDRD@meta.data %>%
  as.data.frame() %>%
  dplyr::arrange(Cell) %>%
  dplyr::group_by(Cell,Genotype) %>% 
  dplyr::summarise(cnt = n()) %>% 
  dplyr::group_by(Cell) %>%  
  dplyr::mutate(percent = 100*round(cnt/sum(cnt), 2)) %>% 
  dplyr::group_by(Cell) %>%
  ggbarplot("Cell", "percent",
            fill = "Genotype", color = "Genotype", palette = color,
            label = TRUE,lab.pos = "in", lab.col = "white") +
  scale_y_continuous(labels = scales::dollar_format(suffix = "%", prefix = "")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="bottom") +
  xlab("")

ggsave("output_NoDRD/06_CellProportion_By_Genotype.pdf", plot = prop_cell, width = 16, height = 18, units = "in", dpi = 150)

# 03 Cell Count
cnt <- table(emilia_noDRD$Cell, emilia_noDRD$Genotype)
openxlsx::write.xlsx(cnt,file = "output_NoDRD/Cell_Count.xlsx", colNames = TRUE, borders = "columns")

#MArkers

all_markers_clustID <- presto::wilcoxauc(emilia_noDRD, 'Cell', assay = 'data')

all_markers.Sign <- all_markers_clustID %>%
  dplyr::filter(padj < 0.05, logFC > 0)

openxlsx::write.xlsx(all_markers.Sign, 
                     file = "output_NoDRD/PrestoByCell_padjLT05_logfcGT0.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")

## top 10 markers for each group
## filter for nominally significant (p<0.05) and over-expressed (auc>0.5)
top10 <- presto::top_markers(all_markers.Sign,n = 10,auc_min = 0.5, pval_max = 0.05)

openxlsx::write.xlsx(top10, 
                     file = "output_NoDRD/PrestoByCell_Top10.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")

emilia_noDRD_slim <- DietSeurat(emilia_noDRD, 
                                         counts = TRUE, 
                                         data = TRUE, 
                                         scale.data = FALSE,
                                         assays="RNA",
                                         dimreducs = c("pca","umap"))

saveRDS(emilia_noDRD_slim, file ="output_NoDRD/06_Final_noDRD_slim.rds")
