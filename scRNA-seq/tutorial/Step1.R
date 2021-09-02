# Macrophagocyte
# initial environment
rm(list=ls())
gc()
setwd("C:/Users/Administrator/Desktop/GSE77288")

# library
library(SingleCellExperiment)
library(ggplot2)
library(scater)

# read data
tung_counts <- read.table("molecules.txt", sep = "\t")
tung_annotation <- read.table("annotation.txt", sep = "\t", header = TRUE)

# note that the data passed to the assay slot has to be a matrix!
tung <- SingleCellExperiment(
  assays = list(counts = as.matrix(tung_counts)),
  colData = tung_annotation
)

# add assay matrix
assay(tung, "logcounts") <- log2(counts(tung) + 1)
# check
logcounts(tung)[1:10, 1:4]

# add annotation
colData(tung)$mean_counts <- colMeans(counts(tung))
colData(tung)$total_counts <- colSums(counts(tung))

# ggplot2 plot
cell_info <- as.data.frame(colData(tung))
ggplot(data = cell_info, aes(x = batch, y = total_counts)) +
  geom_violin(fill = 'brown') + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

# scater plot for specific gene
ggcells(tung, aes(x = batch, y = ENSG00000198938), exprs_values = "logcounts") + 
  geom_violin(fill = 'coral2') + theme_bw() + 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))



