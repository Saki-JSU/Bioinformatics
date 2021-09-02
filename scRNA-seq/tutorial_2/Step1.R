# Step 1:
# initial environment
rm(list=ls())
gc()
setwd("C:/Users/Administrator/Desktop/GSE77288")

# library 
library(scater)
library(data.table)
library(pcaMethods)
library(SC3)
library(scater)
library(SingleCellExperiment)
library(pheatmap)
library(mclust)

# random seed 
set.seed(1234567)

# read data
deng<-readRDS("deng-reads.rds")

# check cell type
table(deng$cell_type2)

# PCA analysis
deng <- runPCA(deng, exprs_values = "logcounts")
dim(reducedDim(deng, "PCA"))
plotPCA(deng, colour_by = "cell_type2")

# estimate value of k
deng <- sc3_estimate_k(deng) 
metadata(deng)$sc3$k_estimation  # number of cluster match with cell_type1
plotPCA(deng, colour_by = "cell_type1")
# SC3 
deng <- sc3(deng, ks = 10, biology = TRUE, n_cores = 1)
# plot consensus matrix
sc3_plot_consensus(deng, k = 10, show_pdata = "cell_type2")
# Silhouette plot
sc3_plot_silhouette(deng, k = 10)
# heatmap plot 
sc3_plot_expression(deng, k = 10, show_pdata = "cell_type2")
# identifier maker genes
sc3_plot_markers(deng, k = 10, show_pdata = "cell_type2")
# PCA with sc3 clusters
plotPCA(deng, colour_by = "sc3_10_clusters")
# compare the results of `SC3` clustering with the original publication cell type labels:
adjustedRandIndex(colData(deng)$cell_type2, colData(deng)$sc3_10_clusters)

# tSNE+k-means clusters
deng <- runTSNE(deng, rand_seed = 1)
plotTSNE(deng)
# k=8 clusters
colData(deng)$tSNE_kmeans <- as.character(kmeans(reducedDims(deng)[[2]], centers = 8)$clust)
plotTSNE(deng, colour_by = "tSNE_kmeans")

# SINCERA
# use the same gene filter as in SC3
input <- logcounts(deng[rowData(deng)$sc3_gene_filter, ])
# perform gene-by-gene per-sample z-score transformation
dat <- apply(input, 1, function(y) scRNA.seq.funcs::z.transform.helper(y))
# hierarchical clustering
dd <- as.dist((1 - cor(t(dat), method = "pearson"))/2)
hc <- hclust(dd, method = "average")
# identifier k
num.singleton <- 0
kk <- 1
for (i in 2:dim(dat)[2]) {
  clusters <- cutree(hc, k = i)
  clustersizes <- as.data.frame(table(clusters))
  singleton.clusters <- which(clustersizes$Freq < 2)
  if (length(singleton.clusters) <= num.singleton) {
    kk <- i
  } else {
    break;
  }
}
cat(kk)
# plot 
pheatmap(
  t(dat),
  cluster_cols = hc,
  cutree_cols = kk,
  kmeans_k = 100,
  show_rownames = FALSE
)

