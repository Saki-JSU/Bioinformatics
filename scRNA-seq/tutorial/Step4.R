# Step 4: Normalisations
# initial environment
rm(list=ls())
gc()
setwd("C:/Users/Administrator/Desktop/GSE77288")

# library
# install.packages("devtools")
# devtools::install_github("hemberg-lab/scRNA.seq.funcs")
library(scRNA.seq.funcs)
library(scater)
library(scran)

# read data
umi <- readRDS("umi.rds")
umi.qc <- umi[! rowData(umi)$discard,! colData(umi)$discard]

# PCA on logcounts_raw data
umi.qc <- runPCA(umi.qc,exprs_values = "logcounts_raw")
plotPCA(umi.qc,colour_by = "batch", size_by = "detected", shape_by = "individual")

# PCA on CPM-normalized data
logcounts(umi.qc) <- log2(calculateCPM(umi.qc) + 1)
umi.qc <- runPCA(umi.qc)
plotPCA(umi.qc, colour_by = "batch", size_by = "detected", shape_by = "individual")
# relative log expression (RLE) plots
plotRLE(umi.qc, exprs_values = "logcounts_raw",colour_by = "batch") + ggtitle("RLE plot for logcounts_raw")
plotRLE(umi.qc, exprs_values = "logcounts",colour_by = "batch") + ggtitle("RLE plot for log2(CPM) counts")

# PCA on scran-normalized Data
# quick cluster
qclust <- quickCluster(umi.qc, min.size = 30)
table(qclust)
# calculate size factor 
umi.qc <- computeSumFactors(umi.qc, clusters = qclust)
umi.qc <- logNormCounts(umi.qc)
# PCA on normalized data
umi.qc <- runPCA(umi.qc)
plotPCA(umi.qc, colour_by = "batch",size_by = "detected", shape_by = "individual")
# RLE plots
plotRLE(umi.qc, exprs_values = "logcounts",colour_by = "batch")
# check to avoid size factor <=0
summary(sizeFactors(umi.qc))

# PCA with Downsampled Data
logcounts(umi.qc) <- log2(Down_Sample_Matrix(counts(umi.qc)) + 1)
umi.qc <- runPCA(umi.qc)
plotPCA(umi.qc,colour_by = "batch",size_by = "detected", shape_by = "individual")
# RLE plots
plotRLE(umi.qc, exprs_values = "logcounts",colour_by = "batch")




















