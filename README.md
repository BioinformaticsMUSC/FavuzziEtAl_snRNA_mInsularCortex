# FavuzziEtAl_snRNA_mInsularCortex
Single Nuclei Mouse Insular Cortex Data Analysis Code

**Sensory neurons gate a cortical plasticity window for immune encoding**
==========================
Activity-dependent plasticity allows brain circuits to encode experience, but this has been established with brief stimuli. Immune responses develop over days and are represented in the insular cortex (IC). How the cortex encodes experiences with this temporal structure is unknown. Here we show that Trpv1+ sensory neurons open a multi-day cortical plasticity window matched to the immune response. Following peripheral immune challenges, Parvalbumin-expressing inhibitory interneurons are activated first, restricting plasticity. Vasoactive intestinal peptide-expressing interneurons respond two days later, disinhibiting IC excitatory neurons. Weaker challenges trigger the first phase but not the second, and preventing the transition blocks immune encoding, whereas lowering inhibition beforehand does not. These results identify a body-to-brain mechanism that gates cortical plasticity on the timescale of slowly evolving physiological experience.
==========================
This repository contains analysis code for the single cell RNA-seq project carried out by researchers 
at the [Favuzzi lab, Yale University](https://www.favuzzilab.com/) and [Berto Lab, MUSC](https://bertolab.org/)

## Cite this

If you use anything in this repository please cite the following publication:

Pre-print URL: https://www.biorxiv.org/content/10.64898/2026.11.02.763631v1

## Access the data with an app:

Here a webapp to analyze the data:

[ShinyApp](https://bioinformatics-musc.shinyapps.io/Favuzzi_snRNA_mm_InsCortex_Preliminary/)

## Files
| directory | contents | code |
| --------- | -------- | -------- |
| [`Initial_Clustering`](Initial_Clustering/) | Code for initial clustering and integration. | 01_Seurat_clustering.R \ 02_Doubletting\ 03_CellCycleScoring.R \ 04_FindMarkers.R |
| [`Enrichment Analysis`](Enrichment/) | Code for Fisher Enrichment | 05_Fisher_EnrichMarkers_AllenPFC_Mouse.R \ 06_Fisher_Enrich_Markers_AllenPFC_Human.R |
| [`Reference Mapping`](RefMapping/) | Code for Azimuth Reference Mapping | 07_Azimuth_ReferenceMapping.R |
| [`Celltype_Annotation`](Final_Labelling/) | Code for Final Labelling | 08_Relable.R |
| [`Classify_Sex`](Classify_Sex/) | Code for sex classification | 09_classifySex.R |
| [`RemoveDRD`](RemoveDRD/) | Code to remove DRD cells and Recluster | 10_RemoveDRD.R |
| [`Celltype_prioritization`](Celltype_prioritization/) | Code for cell type prioritization | 11_AUGUR.R |
| [`DGE`](output_dge/) | Code for DGE analyses | 13_DGE_MaleVsFemale_Rerun.R \ 13_DGE_StimvsCtrl_Rerun.R \ 13a_DGE_Male_DSS_Ctrl.R \ 13b_DGE_Female_DSS_Ctrl.R|
| [`GO`](GO/) | Code for GO | 14_GO_ToppGene.R |
| [`Cell2Cell_Signaling`](Cell2Cell_Signaling/) | Code for cell to cell signaling | 14_GO_ToppGene.R |
