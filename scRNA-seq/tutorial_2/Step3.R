# Step 3: Seurat
# initial environment
rm(list=ls())
gc()
setwd("C:/Users/Administrator/Desktop/GSE77288")

# library 
library(Seurat)
library(ggplot2)
library(SingleR)
library(dplyr)
library(celldex)
library(RColorBrewer)
library(SingleCellExperiment)
library(glmGamPoi)

# read data
adj.matrix <- Read10X("soupX_pbmc10k_filt")
# create Seurate object
srat <- CreateSeuratObject(adj.matrix,project = "pbmc10k") 
srat
# meta data
meta <- srat@meta.data
dim(meta)
head(meta)

# QC 
# calculate MT- percentage
srat[["percent.mt"]] <- PercentageFeatureSet(srat, pattern = "^MT-")
# calculate ribosomal proteins percentage
srat[["percent.rb"]] <- PercentageFeatureSet(srat, pattern = "^RP[SL]")
# doublet annotation
doublets <- read.table("scrublet_calls.tsv",header = F,row.names = 1)
colnames(doublets) <- c("Doublet_score","Is_doublet")
srat <- AddMetaData(srat,doublets)
head(srat[[]])
# violin plots of the selected metadata features
VlnPlot(srat, features = c("nFeature_RNA","nCount_RNA","percent.mt","percent.rb"),ncol = 4,pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))
# plot metadata features against each other and see how they correlate
FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
FeatureScatter(srat, feature1 = "nCount_RNA", feature2 = "percent.rb")
# QC criterion
srat[['QC']] <- ifelse(srat@meta.data$Is_doublet == 'True','Doublet','Pass')
srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < 500 & srat@meta.data$QC == 'Pass','Low_nFeature',srat@meta.data$QC)
srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < 500 & srat@meta.data$QC != 'Pass' & srat@meta.data$QC != 'Low_nFeature',paste('Low_nFeature',srat@meta.data$QC,sep = ','),srat@meta.data$QC)
srat[['QC']] <- ifelse(srat@meta.data$percent.mt > 15 & srat@meta.data$QC == 'Pass','High_MT',srat@meta.data$QC)
srat[['QC']] <- ifelse(srat@meta.data$nFeature_RNA < 500 & srat@meta.data$QC != 'Pass' & srat@meta.data$QC != 'High_MT',paste('High_MT',srat@meta.data$QC,sep = ','),srat@meta.data$QC)
table(srat[['QC']])
# plot gene pass QC
VlnPlot(subset(srat, subset = QC == 'Pass'), 
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.rb"), ncol = 4, pt.size = 0.1) & 
  theme(plot.title = element_text(size=10))

# normalization
srat <- NormalizeData(srat)
# discovers the most variable features (genes)
srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(srat), 10)
top10 
# plot
plot1 <- VariableFeaturePlot(srat)
LabelPoints(plot = plot1, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

# scale
all.genes <- rownames(srat)
srat <- ScaleData(srat, features = all.genes)
# use variable genes to run PCA
srat <- RunPCA(srat, features = VariableFeatures(object = srat))
# visualization
VizDimLoadings(srat, dims = 1:9, reduction = "pca") & 
  theme(axis.text=element_text(size=5), axis.title=element_text(size=8,face="bold"))
DimHeatmap(srat, dims = 1:6, nfeatures = 20, cells = 500, balanced = T)
DimPlot(srat, reduction = "pca")
ElbowPlot(srat)

# cluster
srat <- FindNeighbors(srat, dims = 1:10)
srat <- FindClusters(srat, resolution = 0.5)
# UMP dimension reduction
srat <- RunUMAP(srat, dims = 1:10, verbose = F)
DimPlot(srat,label.size = 4,repel = T,label = T)
# two markers: LILRA4 and TPM2 for DCs, and PPBP and GP1BB for platelets.
FeaturePlot(srat, features = c("LILRA4", "TPM2", "PPBP", "GP1BB"))

# calculate cell cycle scores,
cc.genes.updated.2019
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
srat <- CellCycleScoring(srat, s.features = s.genes, g2m.features = g2m.genes)
table(srat[[]]$Phase)

# SCTransform normalization and clustering
srat <- SCTransform(srat, method = "glmGamPoi", ncells = 8824, 
                    vars.to.regress = c("percent.mt","S.Score","G2M.Score"), verbose = F)
srat
srat <- RunPCA(srat, verbose = F)
srat <- RunUMAP(srat, dims = 1:30, verbose = F)
srat <- FindNeighbors(srat, dims = 1:30, verbose = F)
srat <- FindClusters(srat, verbose = F)
table(srat[[]]$seurat_clusters)