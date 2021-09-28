# pre-process
# initial environment
rm(list=ls())
gc()
setwd("C:/Users/Administrator/Desktop/2021-Pal/code")

# library 
library(Seurat)
library(scater)

# read data
file_name <- list.files("../data/Figure1/")  # file names of all patients
# output
PatientList <- list()

# patient i
for(t in 1:length(file_name)){
  i <- file_name[t]
  # count data
  Counts <- Read10X(paste0("../data/Figure1/",i))  
  # create Seurate object with filter cells and genes
  Seu <- CreateSeuratObject(counts = Counts, min.cells = 0.01 * ncol(Counts))
  Seu <- subset(Seu, subset = nFeature_RNA > 500)
  # filter genes with high mitochondrial genome
  Seu[["percent.mt"]] <- PercentageFeatureSet(Seu, pattern = "^MT-") 
  # filter genes with high nFeature/counts (doublet)
  Seu <- subset(Seu, subset = !isOutlier(Seu$nFeature_RNA) & !isOutlier(Seu$nCount_RNA) & percent.mt < 20)
  # normalization
  Seu <- NormalizeData(Seu)
  # Find Variable Genes
  Seu <- FindVariableFeatures(Seu, selection.method = "vst", nfeatures = 2000)
  # scale data
  Seu <- ScaleData(Seu)
  # Patient ID
  Seu@meta.data$Patient <- i
  # pre/post menopause
  if(t %in% c(1:4, 6, 9:11))
    Seu@meta.data$Menopause <- "Pre"
  else
    Seu@meta.data$Menopause <- "Post"
  # output 
  PatientList[[t]] <- Seu
}

# output
save(PatientList, file = "N-epi.RData")

