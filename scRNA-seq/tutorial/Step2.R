# Step 1: QC
# initial environment
rm(list=ls())
gc()
setwd("C:/Users/Administrator/Desktop/GSE77288")

# library
library(scater)
library(SingleCellExperiment)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(EnsDb.Hsapiens.v86)

# read data
molecules <- read.delim("molecules.txt",row.names=1)
annotation <- read.delim("annotation.txt",stringsAsFactors = T)

# creat SingleCellExperiment object
umi <- SingleCellExperiment(assays = list(counts = as.matrix(molecules)), colData = annotation)
altExp(umi,"ERCC") <- umi[grep("^ERCC-",rownames(umi)), ]  # save ERCC genes solely
umi <- umi[grep("^ERCC-",rownames(umi),invert = T), ]

# map ENSEMBL IDs to gene symbols with org.Hs.eg.db
gene_names <- mapIds(org.Hs.eg.db, keys=rownames(umi), keytype="ENSEMBL", columns="SYMBOL",column="SYMBOL")
rowData(umi)$SYMBOL <- gene_names
# remove gene with no symbols
umi <- umi[! is.na(rowData(umi)$SYMBOL),]
# check mitochondrial proteins
grep("^MT-",rowData(umi)$SYMBOL,value = T) # org.Hs.eg.db does not identify mito, even if they exist
# check ribosomal proteins (which start with RPL or RPS)
grep("^RP[LS]",rowData(umi)$SYMBOL,value = T)

# map ENSEMBL IDs to gene symbols with EnsDb.Hsapiens.v86 (recommanded)
ensdb_genes <- genes(EnsDb.Hsapiens.v86)
MT_names <- ensdb_genes[seqnames(ensdb_genes) == "MT"]$gene_id
is_mito <- rownames(umi) %in% MT_names
table(is_mito)

# add per-cell and per-gene metrics
umi_cell <- perCellQCMetrics(umi,subsets=list(Mito=is_mito))
umi_feature <- perFeatureQCMetrics(umi)
head(umi_cell)




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



