# Step 3
# initial environment
rm(list=ls())
gc()
setwd("C:/Users/Administrator/Desktop/GSE84465")

# library 
library(Seurat)
library(DropletUtils)
library(scDblFinder)
library(scater)
library(celldex)
library(SingleR)

# read data
load("./code/Step2.RData")

# QC for scRNA-seq specifically
# empty droplet check
# e.out <- emptyDrops(GetAssayData(GBM, slot = "counts", assay = "RNA"))  
# it reports an error if this QC has been done
# double droplets check
db.test <- findDoubletClusters(GetAssayData(GBM, slot = "counts", assay = "RNA"), 
                          clusters = GBM@meta.data$seurat_clusters)
chosen.doublet <- rownames(db.test)[isOutlier(db.test$num.de, type = "lower", log = TRUE)]
# empty chosen.doublet means no double droplets

# check cell cycle gene
# example in Seurat 
length(c(cc.genes$s.genes, cc.genes$g2m.genes))
# check cell cycle gene in HVG
CaseMatch(c(cc.genes$s.genes, cc.genes$g2m.genes), VariableFeatures(GBM))
# add into meta data
g2m.genes <- cc.genes$g2m.genes
g2m.genes <- CaseMatch(search = g2m.genes, match = rownames(GBM))
s_genes <- cc.genes$s.genes
s_genes <- CaseMatch(search = s_genes, match = rownames(GBM))
GBM <- CellCycleScoring(object = GBM, g2m.features = g2m.genes, s.features = s_genes)
# check the impact of cell cycle gene in clustering
GBM <- RunPCA(GBM, features = c(s_genes, g2m.genes))
DimPlot(GBM, reduction = "pca", group.by = "Phase")
# if plot are clustering well, it shows the impact of cell cycle gene is minor

# cell annotation by SingleR package
refdata <- HumanPrimaryCellAtlasData()  # reference data
testdata <- GetAssayData(GBM, slot = "data")
clusters <- GBM@meta.data$seurat_clusters
cellpred <- SingleR(test = testdata, ref = refdata, labels = refdata$label.fine,  # fine take longer time
                    clusters = clusters, assay.type.test = "logcounts", assay.type.ref = "logcounts")
# annotation results
cellpred$labels
celltype <- data.frame(ClusterID = rownames(cellpred), celltype = cellpred$labels, stringsAsFactors = F)
# add into Seurat
GBM@meta.data$celltype <- "NA"
for(i in 1:nrow(celltype)){
  GBM@meta.data$celltype[which(GBM@meta.data$seurat_clusters == celltype$ClusterID[i])] <- celltype$celltype[i]
}
# plot
DimPlot(GBM, group.by = "celltype", label = F, reduction = "tsne")
save(GBM, file = "Step3.RData")






