# Step 2: Differential Genes
# initial environment
rm(list=ls())
gc()
setwd("C:/Users/Administrator/Desktop/GSE77288")

# library 
library(scRNA.seq.funcs)
library(edgeR)
library(monocle)
library(MAST)
library(ROCR)

# random seed 
set.seed(1)

# true value
DE <- read.table("TPs.txt")
notDE <- read.table("TNs.txt")
GroundTruth <- list(
  DE = as.character(unlist(DE)), 
  notDE = as.character(unlist(notDE))
)

# read data
molecules <- read.table("molecules.txt", sep = "\t")
anno <- read.table("annotation.txt", sep = "\t", header = TRUE)
# remain comparison data
keep <- anno[,1] == "NA19101" | anno[,1] == "NA19239"  
data <- molecules[,keep]
group <- anno[keep,1]
batch <- anno[keep,4]
# remove genes that aren't expressed in at least 6 cells
gkeep <- rowSums(data > 0) > 5;
counts <- data[gkeep,]
# Library size normalization
lib_size <- colSums(counts)
norm <- t(t(counts)/lib_size * median(lib_size)) 
# Variant of CPM for datasets with library sizes of fewer than 1 mil molecules

# Kolmogorov-Smirnov Test for DE
pVals <- apply(
  norm, 1, function(x) {
    ks.test(
      x[group == "NA19101"], 
      x[group == "NA19239"]
    )$p.value
  }
)
# multiple testing correction
pVals <- p.adjust(pVals, method = "fdr")
# evaluate accuracy
sigDE <- names(pVals)[pVals < 0.05]
length(sigDE) 
# Number of KS-DE genes that are true DE genes
sum(GroundTruth$DE %in% sigDE) 
# Number of KS-DE genes that are truly not-DE
sum(GroundTruth$notDE %in% sigDE)
# TPR and FPR
tp <- sum(GroundTruth$DE %in% sigDE)
fp <- sum(GroundTruth$notDE %in% sigDE)
tn <- sum(GroundTruth$notDE %in% names(pVals)[pVals >= 0.05])
fn <- sum(GroundTruth$DE %in% names(pVals)[pVals >= 0.05])
tpr <- tp/(tp + fn)
fpr <- fp/(fp + tn)
# TPR is much higher than the FPR indicating the KS test is identifying DE genes.
cat(c(tpr, fpr))
# ROC plot
pVals <- pVals[names(pVals) %in% GroundTruth$DE | 
                 names(pVals) %in% GroundTruth$notDE] 
truth <- rep(1, times = length(pVals));
truth[names(pVals) %in% GroundTruth$DE] = 0;
pred <- ROCR::prediction(pVals, truth)
perf <- ROCR::performance(pred, "tpr", "fpr")
ROCR::plot(perf)
# AUC
aucObj <- ROCR::performance(pred, "auc")
aucObj@y.values[[1]] 

# Wilcox/Mann-Whitney-U Test for DE
pVals <- apply(
  norm, 1, function(x) {
    wilcox.test(
      x[group == "NA19101"], 
      x[group == "NA19239"]
    )$p.value
  }
)
# multiple testing correction
pVals <- p.adjust(pVals, method = "fdr")
# AUC
DE_Quality_AUC <- function(pVals, plot=TRUE) {
  pVals <- pVals[names(pVals) %in% GroundTruth$DE | 
                   names(pVals) %in% GroundTruth$notDE]
  truth <- rep(1, times = length(pVals));
  truth[names(pVals) %in% GroundTruth$DE] = 0;
  pred <- ROCR::prediction(pVals, truth)
  perf <- ROCR::performance(pred, "tpr", "fpr")
  if (plot)
    ROCR::plot(perf)
  aucObj <- ROCR::performance(pred, "auc")
  return(aucObj@y.values[[1]])
}
DE_Quality_AUC(pVals)

# edgeR for DE
dge <- DGEList(
  counts = counts, 
  norm.factors = rep(1, length(counts[1,])), 
  group = group
)
group_edgeR <- factor(group)
design <- model.matrix(~ group_edgeR)
dge <- estimateDisp(dge, design = design, trend.method = "none")
fit <- glmFit(dge, design)
res <- glmLRT(fit)
pVals <- res$table[,4]
names(pVals) <- rownames(res$table)
# AUC plot
pVals <- p.adjust(pVals, method = "fdr")
DE_Quality_AUC(pVals)

# Monocle for DE
pd <- data.frame(group = group, batch = batch)
rownames(pd) <- colnames(counts)
pd <- new("AnnotatedDataFrame", data = pd)
Obj <- newCellDataSet(
  as.matrix(counts), 
  phenoData = pd, 
  expressionFamily = negbinomial.size()
)
Obj <- estimateSizeFactors(Obj)
Obj <- estimateDispersions(Obj)
res <- differentialGeneTest(Obj, fullModelFormulaStr = "~group")
pVals <- res[,3]
names(pVals) <- rownames(res)
pVals <- p.adjust(pVals, method = "fdr")
DE_Quality_AUC(pVals)