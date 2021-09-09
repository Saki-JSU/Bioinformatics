# Step 2
# initial environment
rm(list=ls())
gc()
setwd("C:/Users/Administrator/Desktop/GSE84465")

# library 
library(Seurat)
library(ggplot2)
library(tidyverse)

# read data
load("./code/Step1.RData")

# dimension reduction
# normalization
GBM <- ScaleData(GBM, features = rownames(GBM))
# compare 
GetAssayData(GBM, slot = "counts", assay = "RNA")
GetAssayData(GBM, slot = "scale.data", assay = "RNA")
# PCA with HVG
GBM <- RunPCA(GBM, features = VariableFeatures(GBM), seed.use = 3)  # seed.use is used to control the direction of figure
# plot PC1 VS PC2
Fig1D <- DimPlot(GBM, reduction = "pca", group.by = "Patient_ID") + labs(tag = "D")
Fig1D 
# PC selection
GBM <- JackStraw(GBM, reduction = "pca", dims = 20)
GBM <- ScoreJackStraw(GBM, dims = 1:20)
# PCA plot
Fig1E_1 <- JackStrawPlot(GBM, dims = 1:20, reduction = "pca") + 
  theme(legend.position = "bottom") + labs(tag = "E")
Fig1E_1
Fig1E_2 <- ElbowPlot(GBM, ndims = 20, reduction = "pca")
Fig1E_2
# combine plots
Fig1D | (Fig1E_1 | Fig1E_2)

# clustering based on PCA
pc.num <- 1:20
GBM <- FindNeighbors(GBM, dims = pc.num)
GBM <- FindClusters(GBM, resolution = 0.5)
# tSNE dimension reduction
GBM <- RunTSNE(GBM, dims = pc.num)
Fig1F <- DimPlot(GBM, reduction = "tsne", label = TRUE) + labs(tag = "F")
Fig1F

# find marker gene
GBM <- NormalizeData(GBM, normalization.method = "LogNormalize")  # result saved to "data" slot
GBM@assays$RNA@data
diff.wilcox <- FindAllMarkers(GBM)  # default method is wilcox test
# filter marker gene
all.marker <- diff.wilcox %>% select(gene, everything()) %>% 
  subset(p_val<0.05 & abs(diff.wilcox$avg_log2FC) > 0.5)
# output
save(all.marker, file = "markergene.RData")
# each cluster find top 10 marker gene
top10 <- all.marker %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10 <- CaseMatch(search = as.vector(top10$gene), match = rownames(GBM))
length(top10)
length(unique(top10))
# plot
Fig1G <- DoHeatmap(GBM, features = top10, group.by = "seurat_clusters")
Fig1G
Fig1F | Fig1G

# combine figures 
Fig1 <- (Fig1A | Fig1B | Fig1C) /
  (Fig1D | Fig1E) / 
  (Fig1F | Fig1G)

save(GBM, file = "Step2.RData")









