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
  library(speckle)
library(SingleCellExperiment)
library(CellBench)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(cellXY)
})

setwd("/Users/SuganyaSubramanian/Emilia")
emilia_seurat <- readRDS("output_Relabel/02_SeuratObj_30pcs_05res_labelled_withPredictions_slim.rds")
emilia_seurat

count_mxt <- as.matrix(emilia_seurat@assays$RNA@counts)
head(rownames(count_mxt))
head(rownames(emilia_seurat))

# Preprocess data
pro.data <- preprocess(count_mxt, genome="Mm", qc = TRUE)

# Look at counts on Xist and superY.all
pdf("output_Relabel/08_ClassifySex_XistCounts.pdf", width = 8, height = 5)
plot(pro.data$tcm.final$Xist, pro.data$tcm.final$superY)
dev.off()

# Cells that are identified as low quality
dis_cell <- pro.data$discarded.cells
length(dis_cell)

# Cells with zero counts on Xist and all Y genes
zero_cell <- pro.data$zero.cells
length(zero_cell)

sex <- classifySex(count_mxt, genome="Mm") #seurat count matrix is with rownames gene symbols

head(sex$prediction)
head(rownames(sex))

table(sex$prediction)

save(sex,file = "output_Relabel/ClassifySex_meta.RData")

emilia_seurat <- AddMetaData(
  object = emilia_seurat,
  metadata = sex)
head(emilia_seurat)

emilia_seurat$SexPrediction <- emilia_seurat$prediction
emilia_seurat$prediction <- NULL

table(emilia_seurat$SexPrediction)

save(emilia_seurat, file = "output_Relabel/04_Final_slim_with_ClassifySex.RData")

pdf("output_Relabel/09_UMAP_SexPrediction.pdf", width = 8, height = 7)
DimPlot(object = emilia_seurat,group.by = "SexPrediction",
        raster=FALSE, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
dev.off()

#### Prediction based on gene expression

setwd("/Users/SuganyaSubramanian/Emilia")
load("output_Relabel/04_Final_slim_with_ClassifySex.RData")

FeaturePlot(emilia_seurat, feature = c("Eif2s3y","Xist"))

FeaturePlot(emilia_seurat, feature = c("Eif2s3y","Eif2s3x"))

#female specific gene
GENE <- "Xist"
EXPR = GetAssayData(object=emilia_seurat,assay="RNA",slot="data")[GENE,]
emilia_seurat@meta.data$Xist <- as.numeric(EXPR)
max(emilia_seurat@meta.data$Xist)
## Xist no gene expression in our as it lncgene

#male specific gene
GENE <- "Eif2s3y"
EXPR = GetAssayData(object=emilia_seurat,assay="RNA",slot="data")[GENE,]
emilia_seurat@meta.data$Eif2s3y <- as.numeric(EXPR)
max(emilia_seurat@meta.data$Eif2s3y)

#Female Specific Gene
GENE <- "Eif2s3x"
EXPR = GetAssayData(object=emilia_seurat,assay="RNA",slot="data")[GENE,]
emilia_seurat@meta.data$Eif2s3x <- as.numeric(EXPR)
max(emilia_seurat@meta.data$Eif2s3x)


emilia_seurat@meta.data <- emilia_seurat@meta.data %>%
  mutate(Gene_Sex_Prediction = case_when(Eif2s3y > 0 & Eif2s3x == 0 ~ "Male",
                          Eif2s3y == 0 & Eif2s3x > 0 ~ "Female",
                          .default = "NA"))
table(emilia_seurat$Gene_Sex_Prediction)

emilia_seurat@meta.data <- emilia_seurat@meta.data %>%
  mutate(Eif2s3y_Sex_Prediction = case_when(Eif2s3y > 0 & Eif2s3x == 0 ~ "Male",
                          Eif2s3y == 0 ~ "Female",
                          .default = "NA"))

table(emilia_seurat$Eif2s3y_Sex_Prediction)

pdf("output_Relabel/09_UMAP_Eif2s3y_SexPrediction.pdf", width = 8, height = 6)
DimPlot(object = emilia_seurat,group.by = "Eif2s3y_Sex_Prediction",
        raster=FALSE, reduction = "umap", label = FALSE, pt.size = 0.5) #+ theme(legend.position="none")
dev.off()

pdf("output_Relabel/09_UMAP_Eif2s3y_SexPrediction_splitby.pdf", width = 18, height = 6)
DimPlot(object = emilia_seurat,split.by = "Eif2s3y_Sex_Prediction",group.by = "Eif2s3y_Sex_Prediction",
        raster=FALSE, reduction = "umap", label = FALSE, pt.size = 0.5, ncol = 3) #+ theme(legend.position="none")
dev.off()

save(emilia_seurat, file = "output_Relabel/04_Final_slim_with_ClassifySex.RData")
























