# Step 4
# initial environment
rm(list=ls())
gc()
setwd("C:/Users/Administrator/Desktop/GSE84465")

# library 
library(monocle)

# read data
load("./code/Step3.RData")

# pseudo-time analysis
# create CellDataSet object
data <- as(as.matrix(GBM@assays$RNA@counts), "sparseMatrix")
pd <- new("AnnotatedDataFrame", data = GBM@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new("AnnotatedDataFrame", data = fData)
mycds <- newCellDataSet(data, phenoData = pd, featureData = fd, expressionFamily = negbinomial.size())
# pre-process
mycds <- estimateSizeFactors(mycds)
mycds <- estimateDispersions(mycds, cores = 4, relative_expr = TRUE)
# use chosen marker gene
load("./code/markergene.RData")
markers.gene <- all.marker$gene
mycds <- setOrderingFilter(mycds, markers.gene)
# dimension reduction
mycds <- reduceDimension(mycds, max_components = 2, reduction_method = "DDRTree")
mycds <- orderCells(mycds)
# trajectory analysis
p1 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters")
p2 <- plot_cell_trajectory(mycds, color_by = "State")




















