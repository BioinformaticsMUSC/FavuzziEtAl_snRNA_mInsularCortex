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
  library(Seurat)
  library(clustree)
  library(leiden)
  library(data.table)
  library(cowplot)
  library(scDblFinder)
  library(BiocSingular)
  library(RColorBrewer)
  library(MetBrewer)
  library(scales)
  library(ComplexHeatmap)
  library(randomcoloR)
  library(patchwork)
})

load("output_sct_Integrated/09_Final_SeuratObject_WithCellCycle.RData")
table(Idents(seuObject_slim_nodoub))
Idents(seuObject_slim_nodoub) <- "SCT_snn_res.0.5"
seuObject_slim_nodoub@meta.data$seurat_clusters <- seuObject_slim_nodoub@meta.data$SCT_snn_res.0.5


# Create Output Dir
dir.create("output_Relabel")

###Adding Allen Mouse enrichment prediction
#cluster 22 had too many prediction, manually editing to remove duplicates
labels <- read.table("output_sct_Integrated/Labels_Clusters_Res05.txt",header=T,sep="\t")

new <- labels %>% separate(Rows, c("ID", "seurat_clusters"),"_",remove = FALSE)

new <- new[order(as.numeric(as.character(new$seurat_clusters))), ]

current.cluster.ids <- new$seurat_clusters
new.cluster.ids <- as.character(new$variable)

seuObject_slim_nodoub@active.ident <- plyr::mapvalues(x = seuObject_slim_nodoub@active.ident, 
                                                     from = current.cluster.ids, 
                                                     to = new.cluster.ids)

#Add it as ident
seuObject_slim_nodoub@meta.data$Allen_Mouse_Prediction <- seuObject_slim_nodoub@active.ident

#change the ident back to seurat clusters to add human prediction
Idents(seuObject_slim_nodoub) <- "seurat_clusters"

labels <- read.table("output_sct_Integrated/Labels_Clusters_human_Res05.txt",header=T,sep="\t")

new <- labels %>% separate(Rows, c("ID", "seurat_clusters"),"_",remove = FALSE)

new <- new[order(as.numeric(as.character(new$seurat_clusters))), ]

current.cluster.ids <- new$seurat_clusters
new.cluster.ids <- as.character(new$variable)


seuObject_slim_nodoub@active.ident <- plyr::mapvalues(x = seuObject_slim_nodoub@active.ident, 
                                                     from = current.cluster.ids, 
                                                     to = new.cluster.ids)

#Add it as ident
seuObject_slim_nodoub@meta.data$Allen_Human_Prediction <- seuObject_slim_nodoub@active.ident

head(seuObject_slim_nodoub)

#### 3. Add azimuth prediction
predictions <- read.delim('Azimuth/azimuth_pred.tsv', row.names = 1)

seuObject_slim_nodoub <- AddMetaData(
  object = seuObject_slim_nodoub,
  metadata = predictions)

seuObject_slim_nodoub$Azimuth_Prediction <- seuObject_slim_nodoub$predicted.subclass
seuObject_slim_nodoub$predicted.subclass <- NULL

head(seuObject_slim_nodoub)

#change the ident back to seurat clusters to my annotations
Idents(seuObject_slim_nodoub) <- "seurat_clusters"

pdf("output_sct_Integrated/16_UMAP_All_Predictions.pdf", width = 20, height = 6)
p1 <- DimPlot(object = seuObject_slim_nodoub, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seuObject_slim_nodoub, reduction = "umap", label = TRUE, pt.size = 0.5, group.by="Azimuth_Prediction") + theme(legend.position="none")
p3 <- DimPlot(object = seuObject_slim_nodoub, reduction = "umap", label = TRUE, pt.size = 0.5, group.by="Allen_Human_Prediction") + theme(legend.position="none")
p4 <- DimPlot(object = seuObject_slim_nodoub, reduction = "umap", label = TRUE, pt.size = 0.5, group.by="Allen_Mouse_Prediction") + theme(legend.position="none")
plot_grid(p1, p2, p3, p4, ncol = 4)
dev.off()

#### Final Annotation
Idents(seuObject_slim_nodoub) <- "seurat_clusters"

Labels <- read.table("output_sct_Integrated/Final_Label_Clusters.txt",header=T,sep="\t")


#Step2: Get labels
current.cluster.ids <- Labels$seurat_clusters
new.cluster.ids <- as.character(Labels$Cell)

current.cluster.ids
new.cluster.ids


table(Idents(seuObject_slim_nodoub))

##Step 3: Rename cluster with new Labels
seuObject_slim_nodoub@active.ident <- plyr::mapvalues(x = seuObject_slim_nodoub@active.ident, 
                                                           from = current.cluster.ids, 
                                                           to = new.cluster.ids)


### Add New column - Cell with corresponding Labels
seuObject_slim_nodoub@meta.data$Cell <- seuObject_slim_nodoub@active.ident

table(seuObject_slim_nodoub@meta.data$Cell)

seuObject_slim_nodoub@meta.data <- seuObject_slim_nodoub@meta.data %>%
  mutate(Cell_Class = case_when(grepl("Oligo", Cell) ~ "Oligodendrocytes",
                                grepl("OPC", Cell) ~ "OPC",
                                grepl("Astro", Cell) ~ "Astrocytes",
                                grepl("Micro", Cell) ~ "Microglia",
                                #grepl("Neuron", Cell) ~ "Neurons",
                                grepl("IT|CT|STR|ENT|L6b", Cell) ~ "Exc_Neurons",
                                grepl("Inh", Cell) ~ "Inh_Neurons",
                                grepl("Endo", Cell) ~ "Endothelial"))

table(seuObject_slim_nodoub@meta.data$Cell_Class)
Astrocytes      Endothelial      Exc_Neurons      Inh_Neurons        Microglia 
            1246              148            11572             1883              772 
Oligodendrocytes              OPC 
            1111              640 


##Save new Seurat data
save(seuObject_slim_nodoub,file = "output_Relabel/01_SeuratObj_30pcs_05res_labelled_withPredictions.RData")

pdf("output_Relabel/01_UMAP_Labelled.pdf", width = 7, height = 6)
DimPlot(object = seuObject_slim_nodoub, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
dev.off()

pdf("output_Relabel/01_UMAP_Labelled_CellClass.pdf", width = 7, height = 6)
DimPlot(object = seuObject_slim_nodoub, reduction = "umap", label = TRUE, pt.size = 0.5, group.by="Cell_Class") + theme(legend.position="none")
dev.off()

pdf("output_Relabel/02_UMAP_Labelled_Groupby.pdf", width = 14, height = 6)
p1 <- DimPlot(object = seuObject_slim_nodoub, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seuObject_slim_nodoub, reduction = "umap", label = FALSE, pt.size = 0.5, group.by="Genotype")
plot_grid(p1, p2)
dev.off()

pdf("output_Relabel/03_UMAP_Labelled_withclusters.pdf", width = 14, height = 6)
Idents(seuObject_slim_nodoub) <- "seurat_clusters"
p1 <- DimPlot(object = seuObject_slim_nodoub, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
Idents(seuObject_slim_nodoub) <- "Cell"
p2 <- DimPlot(object = seuObject_slim_nodoub, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
plot_grid(p1, p2)
dev.off()

pdf("output_Relabel/04_UMAP_Labelled_Splitby_Genotype.pdf", width = 14, height = 6)
DimPlot(object = seuObject_slim_nodoub, reduction = "umap", label = TRUE,ncol = 2, pt.size = 0.5,split.by = "Genotype") + theme(legend.position="none")
dev.off()

order <- c("Astro_Gja1","Oligo_Mbp","OPC_Vcan","Micro_Dock8",
          "L2/3_IT_Olfm3","L2/3_IT_Lrrtm4","L4/5_IT_Pcsk5","L4/5_IT_Rorb","L4/5_IT_Il1rapl2",
          "L6_CT_Foxp2","L6_CT_Hs3st4","L6b_Ndst4","L6b_ENT_Tshz2","D1/D2_STR",
          "Inh_Pvalb","Inh_Adarb2","Inh_Sst",
          "Endo_Cldn5","VLMC_Slc6a13")

Idents(seuObject_slim_nodoub) <- factor(Idents(seuObject_slim_nodoub), levels = order) # reorder the factor based on cell type

markers <- c("Gja1","Mbp","Vcan","Dock8",
            "Olfm3","Lrrtm4","Pcsk5","Rorb","Il1rapl2",
            "Foxp2","Hs3st4","Ndst4","Tshz2","Drd1","Drd2",
            "Pvalb","Adarb2","Sst",
            "Cldn5","Slc6a13")

pdf("output_Relabel/07_DotPlot_Celltype_Markers.pdf", width = 10, height = 6)
DotPlot(seuObject_slim_nodoub, features = markers) + rotate_x_text(angle = 45)
dev.off()

pdf("output_Relabel/07_DotPlot_seurat_Markers.pdf", width = 10, height = 6)
DotPlot(seuObject_slim_nodoub, features = markers, group.by = "seurat_clusters") + rotate_x_text(angle = 45)
dev.off()

### 2 Cell proportion ####
seuObject_slim_nodoub$Cell <- factor(seuObject_slim_nodoub$Cell, levels = order) # reorder the factor based on cell type

color <- c("#17154f","#b0799a","#e69b00","#355828","#4f2400","#b62684")
#"#e69b00","#355828","#5eb0e5","#ee7762","#4f2400","#b62684")

prop_cell <-  seuObject_slim_nodoub@meta.data %>%
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

ggsave("output_Relabel/06_CellProportion_By_Genotype.pdf", plot = prop_cell, width = 16, height = 18, units = "in", dpi = 150)

# 03 Cell Count
cnt <- table(seuObject_slim_nodoub$Cell, seuObject_slim_nodoub$Genotype)
openxlsx::write.xlsx(cnt,file = "output_Relabel/Cell_Count.xlsx", colNames = TRUE, borders = "columns")

#MArkers

all_markers_clustID <- presto::wilcoxauc(seuObject_slim_nodoub, 'Cell', assay = 'data')

all_markers.Sign <- all_markers_clustID %>%
  dplyr::filter(padj < 0.05, logFC > 0)

openxlsx::write.xlsx(all_markers.Sign, 
                     file = "output_Relabel/PrestoByCell_padjLT05_logfcGT0.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")

## top 10 markers for each group
## filter for nominally significant (p<0.05) and over-expressed (auc>0.5)
top10 <- presto::top_markers(all_markers.Sign,n = 10,auc_min = 0.5, pval_max = 0.05)

openxlsx::write.xlsx(top10, 
                     file = "output_Relabel/PrestoByCell_Top10.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")

seuObject_slim_nodoub_slim <- DietSeurat(seuObject_slim_nodoub, 
                                         counts = TRUE, 
                                         data = TRUE, 
                                         scale.data = FALSE,
                                         assays="RNA",
                                         dimreducs = c("pca","umap"))

saveRDS(seuObject_slim_nodoub_slim, file ="output_Relabel/02_SeuratObj_30pcs_05res_labelled_withPredictions_slim.rds")