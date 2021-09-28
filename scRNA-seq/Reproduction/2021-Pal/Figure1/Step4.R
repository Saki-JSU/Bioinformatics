# integration for multiple patients
# Figure 1I
# initial environment
rm(list=ls())
gc()
setwd("C:/Users/Administrator/Desktop/2021-Pal/code/Figure1/")

# library 
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)

# random seed
set.seed(2021)

# read data
load("N-epi_integrated.RData")
load("N-epi_integrated_cluster.RData")

# clustering
PatientList.integrated <- FindNeighbors(PatientList.integrated, reduction = "pca", dims = 1:10)
PatientList.integrated <- FindClusters(PatientList.integrated, resolution = 0.5)
DimPlot(PatientList.integrated, reduction = "tsne")

# assign cell type 
lineage <- PatientList.integrated@meta.data$seurat_clusters
levels(lineage)[c(1:8,10:11)] <- "black"
levels(lineage)[2] <- "pink"
PatientList.integrated@meta.data$lineage <- lineage
# plot Figure 1I_1
Fig1I_1 <- DimPlot(PatientList.integrated, group.by = "lineage", cols = c("black", "deeppink"), pt.size = 1.5) + NoLegend() + 
  geom_text(x = -25, y = 3, label = "Basal", color = "white", size = 8) + 
  geom_text(x = 10, y = 18, label = "ML", color = "white", size = 8) + 
  geom_text(x = 0, y = -15, label = "LP", color = "white", size = 8)

# plot Figure 1I_2
PatientList.sce$lineage <- lineage
Fig1I_2 <- plotDiffusionMap(PatientList.sce, colour_by = "lineage") + labs(x = "", y = "") + 
  scale_color_manual(values = c("black", "deeppink")) + NoLegend() + 
  geom_text(x = -0.003, y = -0.005, label = "Lineage-primed cells", color = "black", size = 6)

Fig1I_1 + Fig1I_2