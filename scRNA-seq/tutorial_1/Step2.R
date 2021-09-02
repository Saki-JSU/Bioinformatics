# Step 2: QC
# initial environment
rm(list=ls())
gc()
setwd("C:/Users/Administrator/Desktop/GSE77288")

# library
library(scater)
library(scales)
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
head(umi_feature)

# add the metrics to metadata
umi <- addPerCellQC(umi, subsets=list(Mito=is_mito))
umi <- addPerFeatureQC(umi)

# plot
hist(
  umi$total,
  breaks = 100
)
abline(v = 25000, col = "red")

# find low quality genes
# low number of detected genes, but high MT gene percentage, are hallmarks of a low quality cell
reasons <- quickPerCellQC(umi_cell, sub.fields=c("subsets_Mito_percent", "altexps_ERCC_percent"))
colSums(as.matrix(reasons))
umi$discard <- reasons$discard  # add to metadata

# plot
plotColData(umi, x="sum", y="subsets_Mito_percent", colour_by="discard")
plotColData(umi, x="sum", y="detected", colour_by="discard")
plotColData(umi, x="altexps_ERCC_percent", y="subsets_Mito_percent",colour_by="discard")

# plot by batch
plotColData(umi, x="sum", y="detected", colour_by="discard", other_fields = "individual") + 
  facet_wrap(~individual) + scale_x_continuous(labels = unit_format(unit = "k", scale = 1e-3))
plotColData(umi, x="sum", y="detected", colour_by="discard", other_fields = "replicate") + 
  facet_wrap(~replicate)  + scale_x_continuous(labels = unit_format(unit = "k", scale = 1e-3))

# plot highly expressed genes
plotHighestExprs(umi, exprs_values = "counts", 
                 feature_names_to_plot = "SYMBOL", colour_cells_by="detected")


# only remain genes which were detected (expression value > 1) in 2 or more cells
keep_feature <- nexprs(umi,byrow = TRUE,detection_limit = 1) >= 2
rowData(umi)$discard <- ! keep_feature
table(rowData(umi)$discard)

# transform into log2
assay(umi, "logcounts_raw") <- log2(counts(umi) + 1)

# save results
saveRDS(umi, file = "umi.rds")


