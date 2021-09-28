# integration for multiple patients
# Figure 1C & 1D
# initial environment
rm(list=ls())
gc()
setwd("C:/Users/Administrator/Desktop/2021-Pal/code/Figure1/")

# library 
library(Seurat)
library(ggplot2)

# random seed
set.seed(2021)

# read data
load("N-epi.RData")

# select features
features <- SelectIntegrationFeatures(object.list = PatientList)
# scale and PCA for each patient
PatientList <- lapply(X = PatientList, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})

# find anchors
anchors <- FindIntegrationAnchors(object.list = PatientList, reference = c(1, 2), reduction = "rpca",
                                  dims = 1:50)
PatientList.integrated <- IntegrateData(anchorset = anchors, dims = 1:50)
PatientList <- NULL

# integration
PatientList.integrated <- ScaleData(PatientList.integrated, verbose = FALSE)
PatientList.integrated <- RunPCA(PatientList.integrated, verbose = FALSE)
PatientList.integrated <- RunTSNE(PatientList.integrated, dims = 1:50)
# plot Figure 1C
DimPlot(PatientList.integrated, group.by = "Patient") + labs(color = "Samples") + ggtitle("Figure 1C")

# output 
save(PatientList.integrated, file = "N-epi_integrated.RData")

# plot Figure 1D
TSNE.coordinate <- PatientList.integrated@reductions$tsne@cell.embeddings
Pre <- data.frame(TSNE.coordinate[PatientList.integrated$Menopause == "Pre", ])
Post <- data.frame(TSNE.coordinate[PatientList.integrated$Menopause == "Post", ])
Fig1D_1 <- ggplot() + geom_point(data = Pre, aes(x=tSNE_1,y=tSNE_2), col = "blue") + 
  theme_test() + labs(x = "", y = "", title = "Pre-Menopause") + scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  theme(plot.title = element_text(hjust = 0.5))
Fig1D_2 <- ggplot() + geom_point(data = Post, aes(x=tSNE_1,y=tSNE_2), col = "yellow") + 
  theme_test() + labs(x = "", y = "", title = "Post-Menopause") + scale_x_continuous(breaks = NULL) + scale_y_continuous(breaks = NULL) +
  theme(plot.title = element_text(hjust = 0.5))
Fig1D_1 | Fig1D_2


