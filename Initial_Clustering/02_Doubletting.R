suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
  library(emmeans)
  library(SingleCellExperiment)
  library(scater)
  library(BiocParallel)
  library(ggpubr)
  library(Seurat)
  library(clustree)
  library(leiden)
  library(data.table)
  library(cowplot)
  library(scDblFinder)
  library(BiocSingular)
  #library(scds)
})
load("output_sct_Integrated/05_SeuratObj_SCT_30pcs_05res_Slimmed.RData")

source("/Users/SuganyaSubramanian/Utils.R")

#Convert into Single Cell Experiment object
sce <- as.SingleCellExperiment(seuObject_integrated_slim)

# Doubleting by genotype
sce <- scDblFinder(sce,
                   samples="Genotype", 
                   #BPPARAM=MulticoreParam(3),
                   nfeatures = 3000,
                   dims = 30,
                   dbr.sd = 1)

# add in the new column "Doublets" to make if it doublet or singlet
seuObject_integrated_slim@meta.data$Doublets <- sce$scDblFinder.class

pdf("output_sct_Integrated/08_UMAP_SingleDoublet_Distribution.pdf", width = 6, height = 6)
DimPlot(object = seuObject_integrated_slim, reduction = "umap",cols = c("#AED6F1","#CD6155"), label = FALSE, pt.size = 0.5, group.by="Doublets")
dev.off()

save(seuObject_integrated_slim, file = "output_sct_Integrated/05_SeuratObj_SCT_withDoublets.RData")

df <- seuObject_integrated_slim@meta.data %>% as.data.frame()
tmp <- table(df$Doublets,df$Genotype) %>% as.data.frame()

pdf("output_sct_Integrated/09_GeneExp_Singlet-Doublet.pdf")
ggplot(tmp,aes(x=Var2, y=Freq, fill=Var1)) + 
  geom_bar(stat="identity") +
  scale_fill_manual(values=c("#AED6F1","#CD6155")) +
  geom_text(aes(label=Freq), vjust=3.0,size=3.5) +
  ggtitle("Singlet/Doublet by Genotype") +
  xlab("Genotype") +
  ylab("Gene Expression") +
  theme(legend.position="none")
dev.off()

#Subset only singlets (remove Doublets)
seuObject_slim_nodoub <- subset(seuObject_integrated_slim, subset = Doublets == "singlet")

seuObject_slim_nodoub <- processing_seurat_sctransform(seuObject_slim_nodoub, 
                                                       vars_to_regress = c("nCount_RNA","pMito","pRibo"), 
                                                       npcs = 30, 
                                                       res = c(0.3,0.4,0.5,0.6,0.7,0.8))

DefaultAssay(seuObject_slim_nodoub) <- "RNA"
seuObject_slim_nodoub <- NormalizeData(object = seuObject_slim_nodoub, 
                                       normalization.method = "LogNormalize", 
                                       scale.factor = 10000)

save(seuObject_slim_nodoub, file = "output_sct_Integrated/06_SeuratObj_SCT_30pcs_08res_NoDoublet_Reclustered.RData")

pdf("output_sct_Integrated/10_UMAP_NoDoublet_res03.pdf", width = 10, height = 6)
p1 <- DimPlot(object = seuObject_slim_nodoub, group.by = "SCT_snn_res.0.3",reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seuObject_slim_nodoub, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

pdf("output_sct_Integrated/10_UMAP_NoDoublet_res04.pdf", width = 10, height = 6)
p1 <- DimPlot(object = seuObject_slim_nodoub, group.by = "SCT_snn_res.0.4", reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seuObject_slim_nodoub, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

pdf("output_sct_Integrated/10_UMAP_NoDoublet_res05.pdf", width = 10, height = 6)
p1 <- DimPlot(object = seuObject_slim_nodoub, group.by = "SCT_snn_res.0.5",reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seuObject_slim_nodoub, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

pdf("output_sct_Integrated/10_UMAP_NoDoublet_res06.pdf", width = 10, height = 6)
p1 <- DimPlot(object = seuObject_slim_nodoub, group.by = "SCT_snn_res.0.6",reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seuObject_slim_nodoub, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

pdf("output_sct_Integrated/10_UMAP_NoDoublet_res07.pdf", width = 10, height = 6)
p1 <- DimPlot(object = seuObject_slim_nodoub, group.by = "SCT_snn_res.0.7",reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seuObject_slim_nodoub, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

pdf("output_sct_Integrated/10_UMAP_NoDoublet_res08.pdf", width = 10, height = 6)
p1 <- DimPlot(object = seuObject_slim_nodoub, group.by = "SCT_snn_res.0.8",reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seuObject_slim_nodoub, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

#choose Resolution 
#0.5
Idents(seuObject_slim_nodoub) <- "SCT_snn_res.0.5"
seuObject_slim_nodoub@meta.data$seurat_clusters <- seuObject_slim_nodoub@meta.data$SCT_snn_res.0.5

save(seuObject_slim_nodoub, file = "output_sct_Integrated/07_SeuObject_SCT_30pcs_05res_NoDoublet_slim.RData")

pdf("output_sct_Integrated/10_UMAP_NoDoublet_Final.pdf", width = 10, height = 6)
p1 <- DimPlot(object = seuObject_slim_nodoub, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seuObject_slim_nodoub, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

color <- c("#17154f","#b0799a","#e69b00","#355828","#5eb0e5","#ee7762","#4f2400","#b62684")

prop_cell <-  seuObject_slim_nodoub@meta.data %>%
  as.data.frame() %>%
  arrange(seurat_clusters) %>%
  group_by(seurat_clusters,Genotype) %>% 
  summarise(cnt = n()) %>% 
  group_by(seurat_clusters) %>%  
  mutate(percent = 100*round(cnt/sum(cnt), 2)) %>% 
  group_by(seurat_clusters) %>%
  ggbarplot("seurat_clusters", "percent",
            fill = "Genotype", color = "Genotype", palette = color,
            label = TRUE,lab.pos = "in", lab.col = "white") +
  scale_y_continuous(labels = scales::dollar_format(suffix = "%", prefix = "")) +
  rotate_x_text(angle = 45) +
  theme(legend.position="bottom") +
  xlab("")

ggsave("output_sct_Integrated/12_CellProportion_By_Genotype.pdf", plot = prop_cell, width = 25, height = 15, units = "in", dpi = 150)

pdf("output_sct_Integrated/13_UMAP_Nodoublets_SplitBy_Sample.pdf", width = 10, height = 6)
DimPlot(seuObject_slim_nodoub,split.by = "Genotype",ncol = 2) + theme(legend.position="none")
dev.off()