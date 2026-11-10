suppressPackageStartupMessages({
  library(clusterProfiler)  
  library(org.Hs.eg.db)
  library(tidyverse)
  library(DOSE)
  library(ReactomePA)
  library(enrichplot)
  library(ggupset)
  library(GSEABase)
  library(gprofiler2)
  library(msigdbr)
  library(enrichplot)
  library(scToppR)
  library(ggpmisc)
})

setwd("/Users/SuganyaSubramanian/Emilia/")

dge <- readRDS("/Volumes/SS_Backup/Emilia/DGE_NoDRD/DE_DSS_vs_Ctrl.rds")
head(dge)
table(dge$cell_type)

dir.create("GO_ToppGene")

tmp1 <- dge %>%
  group_by(cell_type) %>%
  mutate(Abs = abs(avg_logFC)) %>%
  #dplyr::filter(p_val_adj < 0.05 & Abs > 0.2) %>%
  mutate(Threshold = if_else(p_val_adj < 0.05 & Abs > 0.3, "TRUE","FALSE")) %>%
  mutate(Direction = case_when(avg_logFC > 0.3 & p_val_adj < 0.05 ~ "DSS-UpReg", 
                               avg_logFC < -0.3 & p_val_adj < 0.05 ~ "Ctrl-DownReg"))

table(tmp1$cell_type,tmp1$Direction)

tmp1 <- tmp1 %>% filter(! cell_type %in% c('Endo_Cldn5','VLMC_Slc6a13'))

toppData <-  toppFun(
  tmp1,
  topp_categories = NULL,
  cluster_col = "cell_type",
  gene_col = "gene",
  p_val_col = "p_val_adj", 
  logFC_col = "avg_logFC", 
  min_genes = 5,
  max_genes = 5000,
  max_results = 50
)
head(toppData)
save(toppData, file = "GO_ToppGene/ToppData.RData")

clusters <- unique(tmp1$cell_type)
clusters

### grouped by celltypes not separated by direction
for (i in 1:length(clusters)) {
  print(clusters[i])
  
  p1 <- toppPlot(
    toppData,
    clusters = clusters[i],
    category = 'GeneOntologyBiologicalProcess',
    num_terms = 10,
    p_val_adj = "BH",
    p_val_display = "log",
    y_axis_text_size = 10,
    width = 5,
    height = 6)+ scale_color_gradient(low="blue", high="red")
  p1 <- delete_layers(p1,"GeomSegment")#This removes line segment
  rfn <- gsub("/","-",clusters[i])
  ggsave(p1, filename = paste0('GO_ToppGene/',rfn,'-BP_New.pdf'), height = 8, width = 6)
  
  p2 <- toppPlot(
    toppData,
    clusters = clusters[i],
    category = 'GeneOntologyMolecularFunction',
    num_terms = 10,
    p_val_adj = "BH",
    p_val_display = "log",
    y_axis_text_size = 10,
    save = FALSE,
    width = 5,
    height = 6)+ scale_color_gradient(low="blue", high="red")
  p2 <- delete_layers(p2,"GeomSegment")#This removes line segment
  rfn <- gsub("/","-",clusters[i])
  ggsave(p2, filename = paste0('GO_ToppGene/',rfn,'-MF_New.pdf'), height = 8, width = 6)
}

### if you need to compare upreg to down reg
for (i in 1:length(clusters)) {
  print(clusters[i])
  
  upreg <- tmp1 %>% dplyr::filter(cell_type == clusters[i] & Direction == "DSS-UpReg")
  #print(head(upreg))
  print(length(upreg$gene))
  
  dwnreg <- tmp1 %>% dplyr::filter(cell_type == clusters[i] & Direction == "Ctrl-DownReg")
  #print(head(dwnreg))
  print(length(dwnreg$gene))
  
  updt <-  toppFun(
    upreg,
    topp_categories = NULL,
    cluster_col = "cell_type",
    gene_col = "gene",
    p_val_col = "p_val_adj", 
    logFC_col = "avg_logFC", 
    min_genes = 5,
    max_genes = 5000,
    max_results = 50
  )
  
  dwndt <-  toppFun(
    dwnreg,
    topp_categories = NULL,
    cluster_col = "cell_type",
    gene_col = "gene",
    p_val_col = "p_val_adj",
    logFC_col = "avg_logFC",
    min_genes = 5,
    max_genes = 5000,
    max_results = 50
  )
  
  print("MF")
  p1 <- toppPlot(
    updt,
    category = 'GeneOntologyMolecularFunction',
    num_terms = 10,
    p_val_adj = "BH",
    p_val_display = "log",
    save = FALSE,
    width = 6,
    height = 5
  ) + ggtitle(paste0(clusters[i],"_Upreg"))
  
  p2 <- toppPlot(
    dwndt,
    category = 'GeneOntologyMolecularFunction',
    num_terms = 10,
    p_val_adj = "BH",
    p_val_display = "log",
    save = FALSE,
    #save_dir = "DGE_SubMPs_Pseudobulk/GO_UP",
    width = 6,
    height = 5
  ) + ggtitle(paste0(clusters[i],"_Dwnreg"))
  
  pmrg <- ggarrange(p1, p2, ncol = 2, nrow = 1)
  rfn <- gsub("/","-",clusters[i])
  fnm <- paste0("GO_ToppGene/",rfn,"_Compare_MF.pdf")
  ggsave(pmrg, filename = fnm, height = 6, width = 15)
  
  print("BP")
  p1 <- toppPlot(
    updt,
    category = 'GeneOntologyBiologicalProcess',
    num_terms = 10,
    p_val_adj = "BH",
    p_val_display = "log",
    save = FALSE,
    #save_dir = "DGE_SubMPs_Pseudobulk/GO_UP",
    width = 6,
    height = 5
  )+ ggtitle(paste0(clusters[i],"_Upreg"))
  
  p2 <- toppPlot(
    dwndt,
    category = 'GeneOntologyBiologicalProcess',
    num_terms = 10,
    p_val_adj = "BH",
    p_val_display = "log",
    save = FALSE,
    #save_dir = "DGE_SubMPs_Pseudobulk/GO_UP",
    width = 6,
    height = 5
  )+ ggtitle(paste0(clusters[i],"_Dwnreg"))
  
  pmrg <- ggarrange(p1, p2, ncol = 2, nrow = 1)
  rfn <- gsub("/","-",clusters[i])
  fnm <- paste0("GO_ToppGene/", rfn,"_Compare_BP.pdf")
  ggsave(pmrg, filename = fnm, height = 6, width = 15)
  
  print("CC")
  p1 <- toppPlot(
    updt,
    category = 'GeneOntologyCellularComponent',
    num_terms = 10,
    p_val_adj = "BH",
    p_val_display = "log",
    save = FALSE,
    #save_dir = "DGE_SubMPs_Pseudobulk/GO_UP",
    width = 6,
    height = 5
  )+ ggtitle(paste0(clusters[i],"_Upreg"))
  
  p2 <- toppPlot(
    dwndt,
    category = 'GeneOntologyCellularComponent',
    num_terms = 10,
    p_val_adj = "BH",
    p_val_display = "log",
    save = FALSE,
    #save_dir = "DGE_SubMPs_Pseudobulk/GO_UP",
    width = 6,
    height = 5
  )+ ggtitle(paste0(clusters[i],"_Dwnreg"))
  
  pmrg <- ggarrange(p1, p2, ncol = 2, nrow = 1)
  rfn <- gsub("/","-",clusters[i])
  fnm <- paste0("GO_ToppGene/", rfn,"_Compare_CC.pdf")
  ggsave(pmrg, filename = fnm, height = 6, width = 15)
  
  print("drug")
  p1 <- toppPlot(
    updt,
    category = 'Drug',
    num_terms = 10,
    p_val_adj = "BH",
    p_val_display = "log",
    save = FALSE,
    #save_dir = "DGE_SubMPs_Pseudobulk/GO_UP",
    width = 6,
    height = 5
  )+ ggtitle(paste0(clusters[i],"_Upreg"))
  
  p2 <- toppPlot(
    dwndt,
    category = 'Drug',
    num_terms = 10,
    p_val_adj = "BH",
    p_val_display = "log",
    save = FALSE,
    #save_dir = "DGE_SubMPs_Pseudobulk/GO_UP",
    width = 6,
    height = 5
  )+ ggtitle(paste0(clusters[i],"_Dwnreg"))
  
  pmrg <- ggarrange(p1, p2, ncol = 2, nrow = 1)
  rfn <- gsub("/","-",clusters[i])
  fnm <- paste0("GO_ToppGene/", rfn,"_Compare_Drug.pdf")
  ggsave(pmrg, filename = fnm, height = 6, width = 15)
  
  print('disease')
  p1 <- toppPlot(
    updt,
    category = "Disease",
    num_terms = 10,
    p_val_adj = "BH",
    p_val_display = "log",
    save = FALSE,
    #save_dir = "DGE_SubMPs_Pseudobulk/GO_UP",
    width = 6,
    height = 5
  )+ ggtitle(paste0(clusters[i],"_Upreg"))
  
  p2 <- toppPlot(
    dwndt,
    category = "Disease",
    num_terms = 10,
    p_val_adj = "BH",
    p_val_display = "log",
    save = FALSE,
    #save_dir = "DGE_SubMPs_Pseudobulk/GO_UP",
    width = 6,
    height = 5
  )+ ggtitle(paste0(clusters[i],"_Dwnreg"))
  
  pmrg <- ggarrange(p1, p2, ncol = 2, nrow = 1)
  rfn <- gsub("/","-",clusters[i])
  fnm <- paste0("GO_ToppGene/",rfn,"_Compare_Dis.pdf")
  ggsave(pmrg, filename = fnm, height = 6, width = 15)
}
