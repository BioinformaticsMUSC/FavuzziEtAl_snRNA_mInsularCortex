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
  library(data.table)
  library(cowplot)
  library(scDblFinder)
  library(garnett)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  library(monocle)
})

load("output_sct_Integrated/09_Final_SeuratObject_WithCellCycle_slim.RData")

############ FIND MARKERS #############
all_markers_clustID <- presto::wilcoxauc(seuObject_slim_nodoub_slim, 'seurat_clusters', assay = 'data')
all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")

all_markers.Sign <- all_markers_clustID %>%
  dplyr::filter(padj < 0.05, logFC > 0)

openxlsx::write.xlsx(all_markers.Sign, 
                     file = "output_sct_Integrated/Res05_Presto_Filteredmarkers_padjLT05_logfcGT0.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")

top10 <- presto::top_markers(all_markers.Sign,n = 10,auc_min = 0.5, pval_max = 0.05)

openxlsx::write.xlsx(top10, 
                     file = "output_sct_Integrated/Res05_Top10.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")


all_markers_clustID <- presto::wilcoxauc(seuObject_slim_nodoub_slim, 'SCT_snn_res.0.3', assay = 'data')
all_markers_clustID$group <- paste("Cluster",all_markers_clustID$group,sep="_")

all_markers.Sign <- all_markers_clustID %>%
  dplyr::filter(padj < 0.05, logFC > 0)

openxlsx::write.xlsx(all_markers.Sign, 
                     file = "output_sct_Integrated/Res03_Presto_Filteredmarkers_padjLT05_logfcGT0.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")

## top 10 markers for each group
## filter for nominally significant (p<0.05) and over-expressed (auc>0.5)
top10 <- presto::top_markers(all_markers.Sign,n = 10,auc_min = 0.5, pval_max = 0.05)

openxlsx::write.xlsx(top10, 
                     file = "output_sct_Integrated/Res03_Top10.xlsx", 
                     colNames = TRUE,
                     rowNames = FALSE, 
                     borders = "columns",
                     sheetName="Markers")



## Fisher Enrichment plot
# Differential Expression
load("output_sct_Integrated/09_Final_SeuratObject_WithCellCycle.RData")

Idents(seuObject_slim_nodoub) <- "SCT_snn_res.0.5"
seuObject_slim_nodoub@meta.data$seurat_clusters <- seuObject_slim_nodoub@meta.data$SCT_snn_res.0.5


markers<- presto::wilcoxauc(seuObject_slim_nodoub, 'seurat_clusters', assay = 'data')

markers$group <- paste("Cluster",markers$group,sep="_")

markers_sign <- markers %>%
                dplyr::filter(padj < 0.05, pct_in > 20, pct_out < 20) %>% 
                dplyr::rename(Gene = feature, Cluster = group)

openxlsx::write.xlsx(markers_sign, 
                        file = "output_sct_Integrated/Markers_Significant_Res05.xlsx", 
                        colNames = TRUE, 
                        borders = "columns",
                        overwrite=T)

write.table(markers_sign,"output_sct_Integrated/Markers_Significant_Res05.txt",sep="\t",quote=F,row.names=F)


df <- markers_sign %>% dplyr::select(Gene,Cluster)
write.table(df,"output_sct_Integrated/Markers_Significant_DGE_Res05.txt",sep="\t",quote=F,row.names=F)

df <- read.table("output_sct_Integrated/Markers_Significant_DGE_Res05.txt",header=T,sep="\t")

library(Orthology.eg.db)

M2H <- function(gns){
egs <- mapIds(org.Mm.eg.db, gns, "ENTREZID","SYMBOL")
mapped <- select(Orthology.eg.db, egs, "Homo.sapiens","Mus.musculus")
mapped$HUM <- as.character(mapIds(org.Hs.eg.db, as.character(mapped$Homo.sapiens), "SYMBOL", "ENTREZID"))
mapped$MUS <- gns
mapped
}

res <- M2H(df$Gene)
comb <- merge(df,res,by.x="Gene",by.y="MUS",all=F)

tmp <- comb %>%
       dplyr::select(HUM,Cluster) %>%
       dplyr::rename(Gene = HUM) %>%
       arrange(Cluster)

tmp1 = data.frame(lapply(tmp, as.character), stringsAsFactors=FALSE)
data <- tmp1[!duplicated(tmp1$Gene), ]
data[data == "NULL"] = NA
data <- na.omit(data)

write.table(data,"output_sct_Integrated/Markers_Significant_DGE_HumanID_Res05_nobiomart.txt",sep="\t",quote=F,row.names=F)

human = biomaRt::useMart('ensembl', dataset = 'hsapiens_gene_ensembl')
mouse = biomaRt::useMart('ensembl', dataset = 'mmusculus_gene_ensembl')

MGI = getLDS(attributes = c("mgi_symbol"), 
  filters = "mgi_symbol", values = df$Gene,mart = mouse, 
  attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)


comb <- merge(df,MGI,by.x="Gene",by.y="MGI.symbol",all=F)

tmp <- comb %>%
       dplyr::select(HGNC.symbol,Cluster) %>%
       dplyr::rename(Gene = HGNC.symbol) %>%
       arrange(Cluster)

write.table(tmp,"output_sct_Integrated/Markers_Significant_DGE_HumanID_Res05.txt",sep="\t",quote=F,row.names=F)

markers <- c("Cd4","Sell","Ccr7","Lef1","Tcf7","Klf2","Il7r", #CD4 fmly
             "Cd8a","Gzmk","Cxcr3","Cd69","Nr4a2", #cd8 fmly
             "Prf1","Nkg7","Gzmb","Klrd1","Klrc1",#nk fmly
             "Cd79a","Ms4a1","Cd19",#Bcell
             "Jchain",#"IGHG3","IGHA1","DERL1","FKBP11",#Plasma cell
             "Cd14","S100a8","S100a9", #monocyte Neutrophils
            "C1qa","C1qc","Gpnmb", #macro
             "Itgax","Il3ra",#DC
             "Gata2", "Kit", #Mast cells
             "Gfap","Aqp4", #Astrocytes
              "Vwf","Pecam1", #Endothelial
              "Mef2c","P2ry12", #Microglia
              "Dlx1","Dlx2","Synpr","Cnr1", "Syt1",#Neurons
              "Mog","Mbp", #Oligodendrocytes
              "Acta2","Tagln","Myh1","Myh9",#Mural Cells or Smooth Muscle Cells
              "Olig1","Vcan","Sox5","Pcdh15",#OPC oligodendrocyte percursor cells
              "Slc6a13", "Dcn",#VLMC
              "Sox2","Nes",#NPC neural percursor cells
              "Vim","Col4a1","Cdk4","Nf1","Tme",#mesenchymal stem cells
               "Cldn5", "Ifi27", "Sparcl1","S100b") #endothelial #ependymal

pdf("output_sct_Integrated/18_DotPlot_Immune$Brain.pdf", width = 30, height =15)
DotPlot(seuObject_slim_nodoub_slim, features = unique(markers)) + rotate_x_text(angle = 45)
dev.off()

pdf("output_sct_Integrated/18_DotPlot_Immune$Brain_Res03.pdf", width = 30, height =15)
DotPlot(seuObject_slim_nodoub_slim, features = unique(markers), group.by = "SCT_snn_res.0.3") + rotate_x_text(angle = 45)
dev.off()

gene_list <- c("Mef2c","Slc17a7","Cux2","Calb1","Ntng1","Il1rapl2","Foxp2","Fezf2","Elavl2","Pvalb","Sst","Vip","Cck","Adarb2","Mobp","Cxcl14","C1qc","Flt1","Pdgfra","C1qa","Fn1","Vtn")

pdf("output_sct_Integrated/18_DotPlot_GeneMArkers_Res03.pdf", width = 10, height =8)
DotPlot(seuObject_slim_nodoub_slim, features = unique(gene_list), group.by = "SCT_snn_res.0.3") + rotate_x_text(angle = 45)
dev.off()

pdf("output_sct_Integrated/18_DotPlot_GeneMArkers_Res05.pdf", width = 10, height =8)
DotPlot(seuObject_slim_nodoub_slim, features = unique(gene_list), group.by = "SCT_snn_res.0.5") + rotate_x_text(angle = 45)
dev.off()

gene_list <- c("Gja1","Id3",#Astro
                "Calb1", #Exc_L3
                "Cux2", #Exc_L2
                "Rorb", #Exc_L4? 
                "Il1rapl2", #Exc_L5
                "Foxp2","Fezf2", #Exc_L6
                "Ntng1", #Exc_L4
                "Pvalb","Sst","Lamp5","Vip", #Inh
                "Slc7a11","Cped1", #VLMC
                "Vtn", #Pericytes
                "Flt1","Cldn5", #endothelial
                "Mrc1","C1qc", #Microglia
                "Bcas1",#ImmOligo
                "Mobp", #oligo
                "Vcan") #OPC

pdf("output_sct_Integrated/18_DotPlot_GeneMArkers(2)_Res03.pdf", width = 10, height =8)
DotPlot(seuObject_slim_nodoub_slim, features = unique(gene_list), group.by = "SCT_snn_res.0.3") + rotate_x_text(angle = 45)
dev.off()

pdf("output_sct_Integrated/18_DotPlot_GeneMArkers(2)_Res05.pdf", width = 10, height =8)
DotPlot(seuObject_slim_nodoub_slim, features = unique(gene_list), group.by = "SCT_snn_res.0.5") + rotate_x_text(angle = 45)
dev.off()

allen <- c("Gad1","Gad2","Slc32a1","Slc17a7","Lamp5","Ndnf","Sncg",
                "Vip","Sst","Chodl","Pvalb","Cux2","Rorb","Fezf2",
                "Sulf1","Fam84b","Sla2","Foxp2","Nxph4",
                "Aqp4",
                "Mbp",
                "Cldn5",
                "Ctss","C1qa") #OPC

pdf("output_sct_Integrated/18_DotPlot_GeneMArkers_Allen_Res03.pdf", width = 10, height =8)
DotPlot(seuObject_slim_nodoub_slim, features = unique(allen), group.by = "SCT_snn_res.0.3") + rotate_x_text(angle = 45)
dev.off()

pdf("output_sct_Integrated/18_DotPlot_GeneMArkers_Allen_Res05.pdf", width = 10, height =8)
DotPlot(seuObject_slim_nodoub_slim, features = unique(allen), group.by = "SCT_snn_res.0.5") + rotate_x_text(angle = 45)
dev.off()

load("/Users/SuganyaSubramanian/BrainData/GeneSets_Allen_Mouse_PFC.RData")
head(GeneSets)
hr_df <- data.frame(GeneSets)

sink('/Users/SuganyaSubramanian/BrainData/GeneSets_Allen_Mouse_PFC.txt')
print(GeneSets)
sink()

write.table(as.data.frame(GeneSets),file="/Users/SuganyaSubramanian/BrainData/GeneSets_Allen_Mouse_PFC.csv", quote=F,sep=",",row.names=F)
split_data <- split(GeneSets, f = data$a1)

load("BrainData/GeneSets_Allen_Mouse_PFC.RData")

### for human stuff
inh <- c("SST","NOS1","PENK","GPM6A","MAPK6","FUOM","PKIA", #L1-3 with SST
         "PVALB","GABRA5","TRDN","ISLR2","FILIP1","RGS5","OSTN","UNC5B", #L4-6 with PVALB
         "VIP","INPP5J","SOCS2","HTR2C","SCGN","USP36","ABCA8","RPS15A","CDH3","PRKCG", #L1-3
         "ADARB2","MSMO1", #L1-3
         "PAX6","SYT10","KLHL4","PCDH17","RAMP1","TTN","DBI",#L1-3
         "LAMP5","RELN","LHX6","TRPC3","CHRD") #L1-3 

pdf("output_sct_Integrated/20_DotPlot_AllMarkers_Inh.pdf", width = 10, height = 6)
DotPlot(seu_DURA_withHarmony, features = inh) + rotate_x_text(angle = 45)
dev.off()

exc <- c("FEZF2","SLC17A7","SLC24A2","SEC11C","BTBD11","MKX","HTR2C", "SYT6",   #L4-6
         "HPCAL1","NXPH4",#L4-6
         "THEMIS","CCK", #L4-6
         "NPY", #L1-3
         "CBLN4","GPR83","WFDC2", #L1-3
         "KRT17", #L4-6
         "RORB","VAMP1","CACNA1G","SMARCA2","FAM198B","SMYD2","CD63","CNTNAP5",
         "PHLDB2","OSTN","ARC","TNK2","SNCB") #L4-6
pdf("output_sct_Integrated/20_DotPlot_AllMarkers_EXC.pdf", width = 10, height = 6)
DotPlot(seu_DURA_withHarmony, features = exc) + rotate_x_text(angle = 45)
dev.off()