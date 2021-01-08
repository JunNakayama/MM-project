library(Seurat)
library(SeuratData)
library(SeuratDisk)
library(SeuratWrappers)
library(scater)
library(Matrix)
library(tidyverse)
library(dplyr)
library(tidyr)
library(reshape2)
library(pheatmap)
library(Rtsne)
library(cowplot)
library(ggplot2)
library(ggsci)
library(scales)
library(MAST)
library(DOSE)
library(patchwork)
library(plotly)
library(monocle)
library(MASS)
library(loomR)
library(RColorBrewer)
library(grDevices)
library(colorRamps)
library(data.table)
library(hexbin)





d1 = read.delim("GSE106218_GEO_processed_MM_raw_TPM_matrix.txt")
d2 = read.delim("GSE110499_GEO_processed_MM_raw_TPM_matrix.txt")
d = cbind(d1, d2)

anns = read.csv("metadata.csv")

MM <- CreateSeuratObject(counts = d, project = "MM", annotation = anns)

MM <- NormalizeData(MM, normalization.method = "LogNormalize", scale.factor = 10000)
MM <- FindVariableFeatures(MM)


MM[["percent.mt"]] <- PercentageFeatureSet(MM, pattern = "^mt-")
MM <- ScaleData(MM, vars.to.regress = "percent.mt")
MM <- RunPCA(MM, npcs = 100, ndims.print = 1:5, nfeatures.print = 5)

MM <- FindNeighbors(MM, reduction = "pca", dims = 1:75, nn.eps = 0.5)
MM <- FindClusters(MM, resolution = 3, n.start = 10)


MM <- RunUMAP(MM, dims = 1:75, min.dist = 0.75)

MM@meta.data[,7:10] <- anns[,3:6]

p <- DimPlot(MM, reduction = "umap", pt.size = 0.1, group.by = "orig.ident", label = TRUE) + ggtitle(label = "UMAP")
p <- AugmentPlot(plot = p)
p + NoLegend()


saveRDS(MM, file = "MM.rds")



colpall = c("#FF80FF", "#80FF00", "#0080FF", "#FF0080", "#FFFF00", "#008000", "#FF8080", "#FF0000", "#008080", "#80FF80", "#808000", "#8080FF", "#800000", "#00FF80", "#00FF00", "#8000FF", "#FF00FF", "#808080",
"#000080", "#800080", "#FF8000", "#0000FF", "#000000", "#00FFFF", "#80FFFF", "#FFFF80")




DimPlot(MM, cols = colpall, group.by = "orig.ident", reduction = "umap")
DimPlot(MM, cols = colpall, group.by = "orig.ident", reduction = "umap") + NoLegend()
DimPlot(MM, cols = colpall, group.by = "orig.ident", reduction = "umap", label = TRUE) + NoLegend()
DimPlot(MM, cols = colpall, group.by = "orig.ident", reduction = "umap", pt.size = 2) + NoLegend()


DimPlot(MM, cols = colpall, group.by = "Patient", reduction = "umap")
DimPlot(MM, cols = colpall, group.by = "Patient", reduction = "umap") + NoLegend()
DimPlot(MM, cols = colpall, group.by = "Patient", reduction = "umap", label = TRUE) + NoLegend()

DimPlot(MM, cols = colpall, group.by = "ISS.stage", reduction = "umap")
DimPlot(MM, cols = colpall, group.by = "ISS.stage", reduction = "umap") + NoLegend()
DimPlot(MM, cols = colpall, group.by = "ISS.stage", reduction = "umap", label = TRUE) + NoLegend()

DimPlot(MM, cols = colpall, group.by = "Death", reduction = "umap")
DimPlot(MM, cols = colpall, group.by = "Death", reduction = "umap") + NoLegend()
DimPlot(MM, cols = colpall, group.by = "Death", reduction = "umap", label = TRUE) + NoLegend()



VlnPlot(MM, features = "SORT1", cols = colpall, group.by = "Patient", pt.size = 0)
VlnPlot(MM, features = "SORT1", cols = colpall, group.by = "Patient", pt.size = 0) + NoLegend()
VlnPlot(MM, features = "SORT1", cols = colpall, group.by = "ISS.stage", pt.size = 0)
VlnPlot(MM, features = "SORT1", cols = colpall, group.by = "ISS.stage", pt.size = 0) + NoLegend()

VlnPlot(MM, features = "LAMP2", cols = colpall, group.by = "Patient", pt.size = 0)
VlnPlot(MM, features = "LAMP2", cols = colpall, group.by = "Patient", pt.size = 0) + NoLegend()
VlnPlot(MM, features = "LAMP2", cols = colpall, group.by = "ISS.stage", pt.size = 0)
VlnPlot(MM, features = "LAMP2", cols = colpall, group.by = "ISS.stage", pt.size = 0) + NoLegend()

FeatureScatter(MM, feature1 = "SORT1", feature2 = "LAMP2", cols = colpall, group.by = "Patient")
FeatureScatter(MM, feature1 = "SORT1", feature2 = "LAMP2", cols = colpall, group.by = "ISS.stage")

FeatureScatter(MM, feature1 = "SORT1", feature2 = "LAMP2", cols = colpall, group.by = "ISS.stage")

FeaturePlot(MM, features = c("SORT1", "LAMP2"), cols = c("grey", "red"))
FeaturePlot(MM, features = c("SORT1", "LAMP2"), blend = TRUE)



# module analysis

list <- read_csv("list.csv")
list = unlist(list)
allgene = rownames(MM)
list1 = list(allgene[allgene %in% unlist(list)])

list2 <- read_csv("list2.csv")
list2 = unlist(list2)
allgene = rownames(MM)
list2 = list(allgene[allgene %in% unlist(list2)])

#> list1
#[[1]]
# [1] "ADAM12"   "ADAMTS5"  "CELP"     "CGB7"     "CHL1"     "CLDN22"   "COL6A5"   "CR1"     
# [9] "CYP2A6"   "DCHS2"    "ECE2"     "FBN3"     "FUT6"     "GABRA4"   "GPR26"    "GRIN2B"  
#[17] "KCNIP4"   "KCNK10"   "LAMA3"    "LONRF2"   "LYPD1"    "MUC21"    "MXRA5"    "NBPF14"  
#[25] "NCR3"     "NOS1"     "OCM"      "OR1E1"    "OR2T12"   "PARVA"    "PCDHB8"   "PDE11A"  
#[33] "PHYHIPL"  "PKD2L2"   "PLGLA"    "PRAMEF14" "PTGER3"   "RIMS4"    "SHC3"     "SLC1A2"  
#[41] "SPATA22"  "SPRR2B"   "SPRR2E"   "TMEM155"  "TP53TG5"  "TSSK2"    "USP17L5" 

#> list2
#[[1]]
# [1] "AGT"      "AHNAK"    "ANXA1"    "AP1S2"    "CA2"      "CAV1"     "CD44"     "CD9"     
# [9] "CDH2"     "CLU"      "CST3"     "CTHRC1"   "CYB5A"    "DSG2"     "EYA2"     "FGFR3"   
#[17] "FUCA2"    "GAS6"     "GPAT2"    "GPRC5D"   "HES6"     "IDH2"     "IGF2BP3"  "ITGB7"   
#[25] "ITM2A"    "ITM2C"    "JUP"      "LRIG1"    "MAF"      "MARCKSL1" "MCAM"     "MEST"    
#[33] "PDE3B"    "PDLIM1"   "PEG10"    "PFKP"     "PFN2"     "PLD4"     "PTP4A3"   "PTPRCAP" 
#[41] "RASSF4"   "RCN1"     "SLC16A14" "SULF1"    "TIMP2"    "TMEM173"  "TMSB15A"  "TSPAN7"  
#[49] "TUBB2B"   "TUBB3"   



mod1 <- AddModuleScore(MM, features = list1, pool = allgene, name = "PC1.")
mod1 <- AddModuleScore(mod1, features = list2, pool = allgene, name = "PC2.")


FeaturePlot(mod1, features = "PC1.1", cols = c("grey", "red"))
FeaturePlot(mod1, features = "PC1.1", cols = c("grey", "red2") , pt.size = 2)
FeaturePlot(mod1, features = "PC2.1", cols = c("grey", "red"))
FeaturePlot(mod1, features = "PC2.1", cols = c("grey", "red2") , pt.size = 2)

VlnPlot(mod1, features = "PC1.1", group.by = "orig.ident")
VlnPlot(mod1, features = "PC2.1", group.by = "orig.ident")

VlnPlot(mod1, features = "PC1.1", group.by = "Death")
VlnPlot(mod1, features = "PC2.1", group.by = "Death")

VlnPlot(MM, features = "ITGB7", cols = colpall, group.by = "Death", pt.size = 0) + NoLegend()







