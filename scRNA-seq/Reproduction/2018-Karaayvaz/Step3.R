# Figure 2a & 2b & 2c
# initial environment
rm(list=ls())
gc()
setwd("C:/Users/saki/Desktop/2018-Karaayvaz/code")

# library
library(SingleCellExperiment)
library(ggplot2)
library(patchwork)
library(broman)
library(Rtsne)

# random seed
set.seed(2021)

# load functions
source("funcs.R")
source("funcs_markers.R")

# load results in the last step
load("sceset_ct.RData")
load("anno_colors.RData")
pd_ct <- colData(sceset_ct)
mat_ct <- sceset_ct@assays@data$exprs

## tsne on cell types
# all cells
to_plot_ct <- unique(pd_ct$cell_types_cl_all)
mat_short_ct <- mat_ct[, which(pd_ct$cell_types_cl_all %in% to_plot_ct)]
pd_short_ct <- pd_ct[which(pd_ct$cell_types_cl_all %in% to_plot_ct), ]
tsne_short_ct <- Rtsne(t(mat_short_ct), perplexity = 30)
colnames(tsne_short_ct$Y) <- c("col1", "col2")
tsne_short_ct$Y <- as.data.frame(tsne_short_ct$Y)
tsne_short_ct$Y$cell_types_cl_all <- pd_short_ct$cell_types_cl_all
tsne_short_ct$Y$cell_types_markers <- pd_short_ct$cell_types_markers
tsne_short_ct$Y$patient <- pd_short_ct$patient
# plot Figure 2a
ggplot(tsne_short_ct$Y, aes(x = col1, y = col2, color = factor(cell_types_cl_all, levels = names(anno_colors$tsne)), 
                            shape = patient)) + 
  geom_point(size = 4) + 
  scale_color_manual(values = anno_colors$tsne) +
  labs(col = "patient", x = "Component 1", y = "Component 2", shape = "patient") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# without epithelial
to_plot_ct <- unique(pd_ct$cell_types_cl_all)[-1]
mat_short_ct <- mat_ct[, which(pd_ct$cell_types_cl_all %in% to_plot_ct)]
pd_short_ct <- pd_ct[which(pd_ct$cell_types_cl_all %in% to_plot_ct), ]
tsne_short_ct <- Rtsne(t(mat_short_ct), perplexity = 30)
colnames(tsne_short_ct$Y) <- c("col1", "col2")
tsne_short_ct$Y <- as.data.frame(tsne_short_ct$Y)
tsne_short_ct$Y$cell_types_cl_all <- pd_short_ct$cell_types_cl_all
tsne_short_ct$Y$cell_types_markers <- pd_short_ct$cell_types_markers
tsne_short_ct$Y$patient <- pd_short_ct$patient
# plot Figure 2b
ggplot(tsne_short_ct$Y, aes(x = col1, y = col2, color = factor(cell_types_cl_all, levels = names(anno_colors$tsne)), 
                            shape = patient)) + 
  geom_point(size = 4) + 
  scale_color_manual(values = anno_colors$tsne) +
  labs(col = "patient", x = "Component 1", y = "Component 2", shape = "patient") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))

# epithelial cells
to_plot_ct <- c("epithelial")
mat_short_ct <- mat_ct[, which(pd_ct$cell_types_cl_all %in% to_plot_ct)]
pd_short_ct <- pd_ct[which(pd_ct$cell_types_cl_all %in% to_plot_ct), ]
tsne_short_ct <- Rtsne(t(mat_short_ct), perplexity = 30)
colnames(tsne_short_ct$Y) <- c("col1", "col2")
tsne_short_ct$Y <- as.data.frame(tsne_short_ct$Y)
tsne_short_ct$Y$cell_types_cl_all <- pd_short_ct$cell_types_cl_all
tsne_short_ct$Y$cell_types_markers <- pd_short_ct$cell_types_markers
tsne_short_ct$Y$patient <- pd_short_ct$patient
# Figure 2c
ggplot(tsne_short_ct$Y, aes(x = col1, y = col2, color = factor(patient, levels = names(anno_colors$patient)), 
                            shape = cell_types_cl_all)) + 
  geom_point(size = 4) + 
  scale_color_manual(values = anno_colors$patient) +
  labs(col = "patient", x = "Component 1", y = "Component 2", shape = "cell type")
#dev.off()