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
})
setwd("/Users/SuganyaSubramanian/Emilia")

dir.create("output_CellChat")

#load("output_Relabel/12_SeuObject_SCT_30pcs_05res_NoDoublet_slim_withCellScoring_Final_labelled_Tcell.RData")

emilia_seu <- readRDS("output_NoDRD/06_Final_noDRD_slim.rds")

Ctrl_subset <- subset(emilia_seu, subset = Genotype == "Ctrl")
table(Ctrl_subset$Genotype) # 8283

DSS_subset <- subset(emilia_seu, subset = Genotype == "DSS")
table(DSS_subset$Genotype) #7361

# Fix the levels 
Idents(Ctrl_subset) <- droplevels(Idents(Ctrl_subset))
Ctrl_subset@meta.data$Cell <- droplevels(Ctrl_subset@meta.data$Cell)

Idents(DSS_subset) <- droplevels(Idents(DSS_subset))
DSS_subset@meta.data$Cell <- droplevels(DSS_subset@meta.data$Cell)

# Create cellchat data
Ctrl_subset_cc <- createCellChat(object=Ctrl_subset,group.by = "Cell")
DSS_subset_cc <- createCellChat(object=DSS_subset,group.by = "Cell")

# Add Mouse DB
CellChatDB <- CellChatDB.mouse #CellChatDB.Human
#showDatabaseCategory(CellChatDB)
#Check Database structure
#dplyr::glimpse(CellChatDB$interaction)

CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")

Ctrl_subset_cc@DB <- CellChatDB.use
DSS_subset_cc@DB <- CellChatDB.use

# Function CellChat
processing_cellchat <- function(object) {
  object %>% 
    subsetData() %>% 
    identifyOverExpressedGenes() %>%
    identifyOverExpressedInteractions() %>%
    projectData(PPI.human) %>%  
    computeCommunProb(raw.use = FALSE, population.size = TRUE,nboot=100) %>%
    filterCommunication(min.cells = 20) %>%
    computeCommunProbPathway(thresh = 0.05) %>%
    aggregateNet() %>%
    netAnalysis_computeCentrality(slot.name = "netP")
}

# Run the analysis for each data
Ctrl_subset_cc <- processing_cellchat(Ctrl_subset_cc)
DSS_subset_cc <- processing_cellchat(DSS_subset_cc)

# Save the data
save(Ctrl_subset_cc, DSS_subset_cc, file = "output_CellChat/CellChat_InputData.RData")

list <- list(Ctrl = Ctrl_subset_cc, DSS = DSS_subset_cc)

merge <- mergeCellChat(list, add.names = names(list))

#Compare the total number of interactions and interaction strength
#To answer on question on whether the cell-cell communication is enhanced or 
#not, CellChat compares the the total number of interactions and interaction 
#strength of the inferred cell-cell communication networks from different 
#biological conditions.

pdf("output_CellChat/00_CellChat_Compare_Interaction_Strength.pdf", width = 6, height = 5)
compareInteractions(merge, show.legend = F)#, group = c(1,2)
dev.off()

list[[1]]@netP$pathways
list[[2]]@netP$pathways

pathway <- union(list[[1]]@netP$pathways, list[[2]]@netP$pathways)#, list[[3]]@netP$pathways)
pathway

groupSize <- as.numeric(table(list[[1]]@idents))
par(mfrow = c(1,2), xpd=TRUE)

pdf("output_CellChat/01_netvisual_circle_Ctrl.pdf", width = 15, height = 15)
netVisual_circle(list[[1]]@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
#netVisual_circle(Em_cc@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()

groupSize <- as.numeric(table(list[[2]]@idents))
par(mfrow = c(1,2), xpd=TRUE)

pdf("output_CellChat/01_netvisual_circle_DSS.pdf", width = 15, height = 15)
netVisual_circle(list[[2]]@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
#netVisual_circle(Em_cc@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
dev.off()
#Error in i_set_edge_attr(x, attr(value, "name"), index = value, value = attr(value,  : 
#Oligo_Mbp has all zero - need to check                                                                           Length of new attribute value must be 1 or 178, the number of target edges, not 165

pdf("output_CellChat/02_netvisual_aggregate_Ctrl.pdf", width = 15, height = 15)
netVisual_aggregate(list[[1]], signaling = list[[1]]@netP$pathways, layout = "chord")
dev.off()

pdf("output_CellChat/02_netvisual_aggregate_DSS.pdf", width = 15, height = 15)
netVisual_aggregate(list[[2]], signaling = list[[2]]@netP$pathways, layout = "chord")
dev.off()


ht1 = netAnalysis_signalingRole_heatmap(list[[1]], pattern = "incoming", signaling = pathway, title = names(list)[1], width = 5, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(list[[2]], pattern = "incoming", signaling = pathway, title = names(list)[2], width = 6, height = 10)

pdf("output_CellChat/03_Heatmaps_Incoming.pdf", width = 15, height = 8)
draw(ht1 + ht2,ht_gap = unit(1.5, "cm"), newpage = FALSE)
dev.off()

ht1 = netAnalysis_signalingRole_heatmap(list[[1]], pattern = "outgoing", signaling = pathway, title = names(list)[1], width = 5, height = 10)
ht2 = netAnalysis_signalingRole_heatmap(list[[2]], pattern = "outgoing", signaling = pathway, title = names(list)[2], width = 6, height = 10)

pdf("output_CellChat/03_Heatmaps_Outgoing.pdf", width = 15, height = 8)
draw(ht1 + ht2, ht_gap = unit(1.5, "cm"), newpage = FALSE)
dev.off()

gg = list()
gg[[1]] <- netAnalysis_signalingRole_scatter(list[[1]],title = "ScatterPlot_Ctrl")
gg[[2]] <- netAnalysis_signalingRole_scatter(list[[2]],title = "ScatterPlot_DSS")

pdf("output_CellChat/07_Scatter.pdf", width = 15, height = 6)
patchwork::wrap_plots(plots = gg)
dev.off()

levels(list[[1]]@idents)

#### for loop to print netvisual bubble ####
for(i in 1:nlevels(list[[1]]@idents)){
  #print(i)
  nm = levels(list[[1]]@idents)[i]
  nm = gsub("/","_",nm)
  print(nm)
  
  #print(nm)
  path <- paste0("output_CellChat/Ctrl_",i,"_ligand-receptor_",nm,".pdf",sep="")
  #print(path)
  p1 <- netVisual_bubble(list[[1]], sources.use = i, targets.use = c(1:17), remove.isolate = FALSE)
  #print(p1)
  pdf(file = path, width = 6, height = 4)
  print(p1)
  dev.off()
}

levels(list[[2]]@idents)

#### for loop to print netvisual bubble ####
for(i in 1:nlevels(list[[2]]@idents)){
  #print(i)
  nm = levels(list[[2]]@idents)[i]
  nm = gsub("/","_",nm)
  print(nm)
  path <- paste0("output_CellChat/DSS_",i,"_ligand-receptor_",nm,".pdf",sep="")
  #print(path)
  p1 <- netVisual_bubble(list[[2]], sources.use = i, targets.use = c(1:17), remove.isolate = FALSE)
  #print(p1)
  pdf(file = path, width = 6, height = 4)
  print(p1)
  dev.off()
}

library(Seurat)
library(CellChat)
library(tidyverse)
library(pheatmap)
df <- TSVZ_cc@net$count %>% as.data.frame()
df <- df[,c(1,2)] # select Microglia and MDM
pdf("TSVZ_Microglia_Interaction.pdf",width=4,height=2)
pheatmap::pheatmap(log2(df+1)%>%t(),color=inferno(100),cluster_cols=F,cluster_rows=F)
dev.off()