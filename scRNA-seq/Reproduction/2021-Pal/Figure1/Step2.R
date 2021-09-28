# integration for multiple patients
# Figure 1E & 1F & 1G
# initial environment
rm(list=ls())
gc()
setwd("C:/Users/Administrator/Desktop/2021-Pal/code/Figure1/")

# library 
library(Seurat)
library(ggplot2)
library(scater)

# read data
load("N-epi_integrated.RData")

# clustering
PatientList.integrated <- FindNeighbors(PatientList.integrated, reduction = "pca", dims = 1:10)
PatientList.integrated <- FindClusters(PatientList.integrated, resolution = 0.5)
DimPlot(PatientList.integrated, reduction = "tsne")

# find marker genes
markers <- FindAllMarkers(PatientList.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# label clusters with markers
markers[markers[,7]=="KRT14", ]  # Basal cell 
markers[markers[,7]=="SLPI" | markers[,7]=="PI3", ]  # luminal progenitors (LP)
markers[markers[,7]=="ANKRD30A", ]  # mature luminal (ML)
# assign cell type 
new.cluster.ids <- c("LP", "LP", "Basal", "ML", "ML", "LP",
                     "ML", "LP", "ML", "Other", "Basal")
names(new.cluster.ids) <- levels(PatientList.integrated)
PatientList.integrated <- RenameIdents(PatientList.integrated, new.cluster.ids)
# Figure 1E
cell.num <- table(PatientList.integrated@active.ident)
ClusterLabels = paste("Cluster", names(cell.num), paste0("(n = ", cell.num, ")"))
ClusterBreaks = names(cell.num)
DimPlot(PatientList.integrated, reduction = "tsne", label = TRUE, pt.size = 0.5) + 
  scale_colour_discrete(breaks = ClusterBreaks, labels = ClusterLabels) + 
  scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) + labs(color = "Clusters")

# Figure 1F
library(ggtern)
# clustering with three major clusters
PatientList.integrated$seurat_clusters <- PatientList.integrated@active.ident
# find marker genes
markers.new <- FindAllMarkers(PatientList.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
# classifier
Basal <- colSums(PatientList.integrated@assays$RNA@counts[markers.new[markers.new[,6] == "Basal", 7], ])
LP <- colSums(PatientList.integrated@assays$RNA@counts[markers.new[markers.new[,6] == "LP", 7], ])
ML <- colSums(PatientList.integrated@assays$RNA@counts[markers.new[markers.new[,6] == "ML", 7], ])
# ternary plot for three major clusters
TernaryData <- data.frame(Basal, LP, ML)
TernaryData <- TernaryData / rowSums(TernaryData)
Clusters <- PatientList.integrated$seurat_clusters
TernaryData <- cbind(TernaryData, Clusters)
TernaryData <- TernaryData[TernaryData[,4] != "Other", ]
ggtern(TernaryData, aes(x = LP, y = Basal, z = ML)) +  geom_point(aes(colour = Clusters)) + NoLegend() + 
  theme_bw() + theme_hidegrid() + theme_hidelabels() + NoLegend()

#  Seurat to SCE
PatientList.sce <- as.SingleCellExperiment(PatientList.integrated)
PatientList.integrated <- NULL
# diffusionmap
PatientList.sce <- runDiffusionMap(PatientList.sce)  # take long time
# plot Figure 1G
plotDiffusionMap(PatientList.sce, colour_by = "ident") + labs(x = "", y = "") + 
  geom_text(x = -0.015, y = -0.0045, label = "Basal") + 
  geom_text(x = 0.005, y = -0.01, label = "ML") + 
  geom_text(x = 0.004, y = 0.005, label = "LP")

# output
save(PatientList.sce, file = "N-epi_integrated_cluster.RData")


