# Step 1
# initial environment
rm(list=ls())
gc()
setwd("C:/Users/Administrator/Desktop/GSE84465")

# library 
library(Seurat)
library(ggplot2)

# read data
raw.data <- read.csv("./data/GSE84465_GBM_All_data.csv", sep = " ")
# delete last 5 rows which are not gene id
# these missing come from htseq software
raw.data <- raw.data[1:(nrow(raw.data)-5),]

# load meta data
meta.data <- read.table("./data/SraRunTable.txt", sep = ",", header = TRUE)
meta.data <- meta.data[, c("plate_id","Well","Tissue","patient_id")]
meta.GBM <- meta.data[meta.data$Tissue=="Tumor",]
GBM_sample <- paste0("X", meta.GBM$plate_id, ".",meta.GBM$Well)
GBM.data <- raw.data[,GBM_sample]

# create meta information
GBM.meta <- data.frame(Patient_ID=meta.GBM$patient_id, row.names = GBM_sample)
# create seurat object with QC
GBM <- CreateSeuratObject(counts = GBM.data, meta.data = GBM.meta, min.cells = 3, min.features = 50)
# filter mitochondrial gene
GBM[["percent.mt"]] <- PercentageFeatureSet(GBM, pattern = "^MT-")  # no mito
pctMT <- 5
GBM <- subset(GBM, subset = percent.mt < pctMT) 
# filter ERCC gene
GBM[["percent.ERCC"]] <- PercentageFeatureSet(GBM, pattern = "^ERCC-")  
pctERCC <- 40  # this threshold should be about 10
GBM <- subset(GBM, subset = percent.ERCC < pctERCC) 

# visualization
# Fig1A Violin Plot
col.num <- length(unique(GBM@meta.data$Patient_ID))
Fig1A.1 <- VlnPlot(GBM, features = c("nFeature_RNA"), group.by = "Patient_ID", cols = rainbow(col.num), pt.size = 2) + 
  theme(legend.position = "none") + labs(tag = "A")
Fig1A.2 <- VlnPlot(GBM, features = c("nCount_RNA"), group.by = "Patient_ID", cols = rainbow(col.num), pt.size = 2) + 
  theme(legend.position = "none")
Fig1A <- Fig1A.1 | Fig1A.2
Fig1A
# Fig1B Scatter Plot
# title is Pearson's correlation of nCount_RNA and nFeature_RNA
Fig1B <- FeatureScatter(GBM, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", group.by = "Patient_ID", pt.size = 1.3) +
  labs(tag = "B")
# Fig1C Find Highly Variable Gene (HVG)
GBM <- FindVariableFeatures(GBM, selection.method = "vst", nfeatures = 1500)
top10 <- head(VariableFeatures(GBM), 10)
plot1 <- VariableFeaturePlot(GBM)
Fig1C <- LabelPoints(plot = plot1, points = top10, repel = TRUE, size = 2.5) +
  theme(legend.position = c(0.1,0.8)) + labs(tag = "C")
Fig1C
# combine Fig1A, Fig1B, Fig1C
Fig1A | Fig1B | Fig1C
# output
save(GBM, file = "Step1.RData")

# following codes are unnecessary for this paper, but may be used in elsewhere
# transform into SCE format
library(scater)
counts <- data.frame(GBM@assays$RNA@counts)
GBM.SCE <- SingleCellExperiment(assays = list(counts = GBM@assays$RNA@counts),
                                colData = GBM@meta.data)
# normalization
stand_exprs(GBM.SCE) <- log2(calculateCPM(GBM.SCE)+1)  # CPM normalization
GBM.SCE <- logNormCounts(GBM.SCE)  # log counts normalization
# plot top10 in four samples
plotExpression(GBM.SCE, top10, x = "Patient_ID", colour_by = "Patient_ID", exprs_values = "logcounts")
p1 <- plotHighestExprs(GBM.SCE, exprs_values = "logcounts")  # default option is top 50
ggsave("./HighestExpr.pdf", plot = p1, width = 15, height = 18)  # suggested to save as a pdf file


