library(RCurl)
library(AnnotationHub)
library(ensembldb)

load("output_sct_Integrated/06_SeuratObj_SCT_30pcs_08res_NoDoublet_Reclustered.RData")

load("/Users/SuganyaSubramanian/MouseCellCycleGenes.rda")

head(g2m_genes)
head(s_genes)

seuObject_slim_nodoub<- CellCycleScoring(seuObject_slim_nodoub,
                                               g2m.features = g2m_genes,
                                               s.features = s_genes)

save(seuObject_slim_nodoub, file ="output_sct_Integrated/09_Final_SeuratObject_WithCellCycle.RData")

seuObject_slim_nodoub_slim <- DietSeurat(seuObject_slim_nodoub, 
                                         counts = TRUE, 
                                         data = TRUE, 
                                         scale.data = FALSE,
                                         assays="RNA",
                                         dimreducs = c("pca","umap"))

save(seuObject_slim_nodoub_slim, file ="output_sct_Integrated/09_Final_SeuratObject_WithCellCycle_slim.RData")

pdf("output_sct_Integrated/14_UMAP_Cellscore.pdf", width = 6, height = 5)
DimPlot(seuObject_slim_nodoub_slim,
        reduction = "umap",
        group.by= "Phase")
dev.off()