# integration for multiple patients
# Figure 1H
# initial environment
rm(list=ls())
gc()
setwd("C:/Users/Administrator/Desktop/2021-Pal/code/Figure1/")

# library 
library(Seurat)
library(ggplot2)
library(cowplot)
library(Matrix)

# read data
load("./N-epi_integrated.RData")

# Figure 1H
PatientList.integrated@assays$integrated@data <- PatientList.integrated@assays$integrated@counts # release memory
PatientList.integrated@assays$RNA@data <- PatientList.integrated@assays$integrated@counts  # release memory
PatientList.integrated@assays$integrated@data <- Matrix(log10(PatientList.integrated@assays$RNA@counts + 1), sparse = TRUE) 

# select marker genes
MarkerGene <- c("EPCAM", "KRT5", "FOXA1", "KIT", "SOX10")
# feature plot for each gene
p <- FeaturePlot(PatientList.integrated, features = MarkerGene, cols = c("grey","red"), keep.scale = "all", combine = FALSE)
# extract the legend from one of the plots
EPCAM <- PatientList.integrated@assays$integrated@data["EPCAM", ]
p.temp <- p[[1]] + theme(legend.direction = "horizontal") + 
  guides(color = guide_colorbar(title = "Log expression", title.position = "top")) + 
  scale_color_gradient2(breaks = c(0, max(EPCAM)), labels = c(0,'Highest'), mid = "grey", high = "red")
legend <- get_legend(p.temp)
# delete legend for each subfigure
for(i in 1:length(p)) {
  p[[i]] <- p[[i]] + NoLegend() + NoAxes()
}
# plot Figure 1H
plot_grid(p[[1]], p[[2]], p[[3]], p[[4]], p[[5]], legend, ncol = 3) 
