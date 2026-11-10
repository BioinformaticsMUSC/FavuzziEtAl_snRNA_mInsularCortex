
predictions <- read.delim('Azimuth/azimuth_pred.tsv', row.names = 1)

seuObject_slim_nodoub_slim <- AddMetaData(
  object = seuObject_slim_nodoub_slim,
  metadata = predictions)

projected.umap <- readRDS('Azimuth/azimuth_umap.Rds')
seuObject_slim_nodoub_slim <- seuObject_slim_nodoub_slim[, Cells(projected.umap)]
seuObject_slim_nodoub_slim[['umap.proj']] <- projected.umap

pdf("output_sct_Integrated/15_UMAP_Azimuth_Prediction.pdf", width = 10, height = 6)
p1 <- DimPlot(object = seuObject_slim_nodoub_slim, reduction = "umap", label = TRUE, pt.size = 0.5) + theme(legend.position="none")
p2 <- DimPlot(object = seuObject_slim_nodoub_slim, reduction = "umap", label = TRUE, pt.size = 0.5, group.by="predicted.subclass") + theme(legend.position="none")
plot_grid(p1, p2)
dev.off()

