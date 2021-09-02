# Step 3: dimension reduction
# initial environment
rm(list=ls())
gc()
setwd("C:/Users/Administrator/Desktop/GSE77288")

# library
library(scater)
library(data.table)

# read data
umi <- readRDS("umi.rds")
umi.qc <- umi[! rowData(umi)$discard,! colData(umi)$discard]

# PCA before QC
# count
umi <- runPCA(umi, exprs_values = "counts")
dim(reducedDim(umi, "PCA"))
plotPCA(umi, colour_by = "batch", size_by = "detected", shape_by = "individual")
# logcount
umi <- runPCA(umi, exprs_values = "logcounts_raw")
dim(reducedDim(umi, "PCA"))
plotPCA(umi, colour_by = "batch", size_by = "detected", shape_by = "individual")
# PCA after QC
umi.qc <- runPCA(umi.qc, exprs_values = "logcounts_raw")
dim(reducedDim(umi.qc, "PCA"))
plotPCA(umi.qc, colour_by = "batch", size_by = "detected", shape_by = "individual")

# TSNE before QC
set.seed(123456)
umi <- runTSNE(umi, exprs_values = "logcounts_raw", perplexity = 130)
plotTSNE(umi, colour_by = "batch", size_by = "detected", shape_by = "individual")
# after QC
set.seed(123456)
umi.qc <- runTSNE(umi.qc, exprs_values = "logcounts_raw", perplexity = 130)
plotTSNE(umi.qc, colour_by = "batch", size_by = "detected", shape_by = "individual")
 
# Detected genes
logcounts(umi.qc) <- assay(umi.qc, "logcounts_raw")
getExplanatoryPCs(umi.qc,variables = "sum")
plotExplanatoryPCs(umi.qc,variables = "sum") 
logcounts(umi.qc) <- NULL

# Explanatory Variables
plotExplanatoryVariables(umi.qc,exprs_values = "logcounts_raw",
                         variables = c("detected","sum","batch",
                                       "individual","altexps_ERCC_percent","subsets_Mito_percent"))



