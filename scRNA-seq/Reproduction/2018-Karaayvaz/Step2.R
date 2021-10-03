# Figure 1c & 1d
# initial environment
rm(list=ls())
gc()
setwd("C:/Users/saki/Desktop/2018-Karaayvaz/code")

# library
library(SingleCellExperiment)
library(broman)
library(ggplot2)
library(tibble)
library(monocle)
library(patchwork)

# load functions
source("funcs.R")
source("funcs_markers.R")

# load results in the last step
load("markers.RData")
load("sceset_final.RData")
load("anno_colors.RData")
mat_norm <- sceset_final@assays@data$exprs
pd_norm <- sceset_final@colData
patients_now <- unique(sceset_final$patient)

## cell types by markers (done as described in SI of paper)
thresh <- 1
cells_markers <- lists_markers(mat_norm, thresh, markers)
epithelial_markers <- cells_markers$epithelial_cells
is_epithelial <- decide_is_epithelial(epithelial_markers)
immune_markers <- cells_markers$immune_cells
is_immune <- decide_is_immune(immune_markers)
other_markers <- cells_markers$other_cells
is_other <- decide_is_other(other_markers)

one_epithelial_marker <- expression_one_epithelial_marker(mat_norm, pd_norm, is_epithelial, epithelial_markers, "pats", 0.5)
is_epithelial[which(one_epithelial_marker$is_epithelial_extra == 1)] <- 1

is_epithelial_simple <- is_epithelial
is_epithelial_simple[which(is_epithelial == 1)] <- "epithelial"
is_immune_simple <- is_immune
is_immune_simple[which(is_immune == "immune_mix")] <- 0
is_other_simple <- is_other
is_other_simple[which(is_other == "other_mix")] <- 0

cells_types <- paste(is_epithelial_simple, is_immune_simple, is_other_simple, sep = "_")
names(cells_types) <- names(is_epithelial)
cell_types <- sapply(strsplit(cells_types, "_"), function(x){
  # none of the cell types (epithelial, immune, other)
  if (sum(x == 0) == 3) return("unknown") else 
    if (sum(x == 0) == 2) return(setdiff(x, "0")) else
      if (sum(c("epithelial", "stroma", "0") %in% x) == 3) return("epithelial") else
        return(paste(setdiff(x, "0"),collapse = "_"))})
cell_types_simple <- cell_types
cell_types_simple[which(sapply(strsplit(cell_types, "_"), length) > 1)] <- "undecided"
table(cell_types_simple)
# update colData and pd_norm
colData(sceset_final)$cell_types_markers <- cell_types_simple
pd_norm <- colData(sceset_final)

## cell types by unsupervised clustering (done as described in SI of paper)
HSMM_clustering_ct <- monocle_unsup_clust_plots(sceset_obj = sceset_final, mat_to_cluster = mat_norm, anno_colors = anno_colors, 
                                                name_in_phenodata = "cluster_allregr_disp", disp_extra = 1, save_plots = 0,
                                                path_plots = NULL, type_pats = "allpats", regress_pat = 1)
# HSMM_clustering_ct$Cluster <- HSMM_clustering_ct$cluster_allregr_disp
# table(HSMM_clustering_ct$Cluster)
# due to changes in Monocle's functions (reduceDimension and clusterCells), the resulting clustering is slightly different from
# the original clustering from the paper. for reproducibility, we read in the original cluster assignment
original_clustering <- readRDS(file = "../data/original_clustering.RDS")
HSMM_clustering_ct$Cluster <- original_clustering
table(HSMM_clustering_ct$Cluster)
# update the sceset_final object
colData(sceset_final) <- cbind(colData(sceset_final), pData(HSMM_clustering_ct)[,c(96:104)])
colData(sceset_final)$cell_types_cl_all <- colData(sceset_final)$cell_types_markers
pd_norm <- colData(sceset_final)

## assign unknown and undecided cell types from clusters
cells_to_assign <- list()
cells_to_reassign <- list()
mean_exprs <- list()
mean_reassign_exprs <- list()
# cluster1: only epithelial cells, so nothing to be done
cluster_here <- 1
table(pd_norm$cell_types_markers[which(pd_norm$Cluster == cluster_here)])
mean_exprs[[cluster_here]] <- NULL
mean_reassign_exprs[[cluster_here]] <- NULL
cells_to_assign[[cluster_here]] <- NULL
cells_to_reassign[[cluster_here]] <- NULL
# cluster2: try to assign the unknown and undecided cells to macrophages
cluster_here <- 2
table(pd_norm$cell_types_markers[which(pd_norm$Cluster == cluster_here)])
to_assign_cluster <- c("undecided", "unknown")
mean_exprs[[cluster_here]] <- compute_mean_expr_types(types = c("epithelial", "immune", "other"), mat_expr = mat_norm, 
                                                      cells_pos = intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% to_assign_cluster)),
                                                      epithelial_markers = epithelial_markers, immune_markers = immune_markers, other_markers = other_markers)
# only assign the cells whose mean immune expression is highest
cells_to_assign[[cluster_here]] <- intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% to_assign_cluster))[which(apply(mean_exprs[[cluster_here]], 1, function(x){if (x[2] >= max(x)) return(1) else {return(0)}}) == 1)]
pd_norm$cell_types_cl_all[cells_to_assign[[cluster_here]]] <- "macrophage"
table(pd_norm$cell_types_cl_all[which(pd_norm$Cluster == cluster_here)])
# also check whether the epithelial cells have high immune expression
mean_reassign_exprs[[cluster_here]] <- compute_mean_expr_types(types = c("epithelial", "immune", "other"), mat_expr = mat_norm, 
                                                               cells_pos = intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% "epithelial")),
                                                               epithelial_markers = epithelial_markers, immune_markers = immune_markers, other_markers = other_markers)
# only re-assign the epithelial cells whose mean immune expression is highest
cells_to_reassign[[cluster_here]] <- intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% "epithelial"))[which(apply(mean_reassign_exprs[[cluster_here]], 1, function(x){if (x[2] >= max(x)) return(1) else {return(0)}}) == 1)]
pd_norm$cell_types_cl_all[cells_to_reassign[[cluster_here]]] <- "macrophage"
table(pd_norm$cell_types_cl_all[which(pd_norm$Cluster == cluster_here)])
# cluster3: mixed, so nothing to be done
cluster_here <- 3
table(pd_norm$cell_types_markers[which(pd_norm$Cluster == cluster_here)])
mean_exprs[[cluster_here]] <- NULL
mean_reassign_exprs[[cluster_here]] <- NULL
cells_to_assign[[cluster_here]] <- NULL
cells_to_reassign[[cluster_here]] <- NULL
# cluster 4: try to assign the unknown and undecided cells to epithelial
cluster_here <- 4
table(pd_norm$cell_types_markers[which(pd_norm$Cluster == cluster_here)])
to_assign_cluster <- c("undecided", "unknown")
mean_exprs[[cluster_here]] <- compute_mean_expr_types(types = c("epithelial", "immune", "other"), mat_expr = mat_norm, 
                                                      cells_pos = intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% to_assign_cluster)),
                                                      epithelial_markers = epithelial_markers, immune_markers = immune_markers, other_markers = other_markers)
# only assign the cells whose epithelial expression is highest
cells_to_assign[[cluster_here]] <- intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% to_assign_cluster))[which(apply(mean_exprs[[cluster_here]], 1, function(x){if (x[1] >= max(x)) return(1) else {return(0)}}) == 1)]
pd_norm$cell_types_cl_all[cells_to_assign[[cluster_here]]] <- "epithelial"
table(pd_norm$cell_types_cl_all[which(pd_norm$Cluster == cluster_here)])
# also check whether the stroma cells have higher epithelial expression
mean_reassign_exprs[[cluster_here]] <- compute_mean_expr_types(types = c("epithelial", "immune", "other"), mat_expr = mat_norm, 
                                                               cells_pos = intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% "stroma")),
                                                               epithelial_markers = epithelial_markers, immune_markers = immune_markers, other_markers = other_markers)
# only re-assign the stroma cells if their epithelial expression is higher
cells_to_reassign[[cluster_here]] <- intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% "stroma"))[which(apply(mean_reassign_exprs[[cluster_here]], 1, function(x){if (x[1] >= max(x)) return(1) else {return(0)}}) == 1)]
pd_norm$cell_types_cl_all[cells_to_reassign[[cluster_here]]] <- "epithelial"
table(pd_norm$cell_types_cl_all[which(pd_norm$Cluster == cluster_here)])
# cluster 5: mixed, so nothing to be done
cluster_here <- 5
table(pd_norm$cell_types_markers[which(pd_norm$Cluster == cluster_here)])
mean_exprs[[cluster_here]] <- NULL
mean_reassign_exprs[[cluster_here]] <- NULL
cells_to_assign[[cluster_here]] <- NULL
cells_to_reassign[[cluster_here]] <- NULL
# cluster 6: try to assign the unknown cells to epithelial
cluster_here <- 6
table(pd_norm$cell_types_markers[which(pd_norm$Cluster == cluster_here)])
to_assign_cluster <- c("unknown")
mean_exprs[[cluster_here]] <- compute_mean_expr_types(types = c("epithelial", "immune", "other"), mat_expr = mat_norm, 
                                                      cells_pos = intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% to_assign_cluster)),
                                                      epithelial_markers = epithelial_markers, immune_markers = immune_markers, other_markers = other_markers)
# only assign the cells whose epithelial expression is highest
cells_to_assign[[cluster_here]] <- intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% to_assign_cluster))[which(apply(mean_exprs[[cluster_here]], 1, function(x){if (x[1] >= max(x)) return(1) else {return(0)}}) == 1)]
pd_norm$cell_types_cl_all[cells_to_assign[[cluster_here]]] <- "epithelial"
table(pd_norm$cell_types_cl_all[which(pd_norm$Cluster == cluster_here)])
# also check whether the Bcell has higher epithelial expression
mean_reassign_exprs[[cluster_here]] <- compute_mean_expr_types(types = c("epithelial", "immune", "other"), mat_expr = mat_norm, 
                                                               cells_pos = intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% "Bcell")),
                                                               epithelial_markers = epithelial_markers, immune_markers = immune_markers, other_markers = other_markers)
# only re-assign the Bcell if its epithelial expression is higher
cells_to_reassign[[cluster_here]] <- intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% "Bcell"))[which(apply(mean_reassign_exprs[[cluster_here]], 1, function(x){if (x[1] >= max(x)) return(1) else {return(0)}}) == 1)]
pd_norm$cell_types_cl_all[cells_to_reassign[[cluster_here]]] <- "epithelial"
table(pd_norm$cell_types_cl_all[which(pd_norm$Cluster == cluster_here)])
# cluster 7: try to assign the unknown cell to epithelial
cluster_here <- 7
table(pd_norm$cell_types_markers[which(pd_norm$Cluster == cluster_here)])
to_assign_cluster <- c("unknown")
mean_exprs[[cluster_here]] <- compute_mean_expr_types(types = c("epithelial", "immune", "other"), mat_expr = mat_norm, 
                                                      cells_pos = intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% to_assign_cluster)),
                                                      epithelial_markers = epithelial_markers, immune_markers = immune_markers, other_markers = other_markers)
# only assign the cells whose epithelial expression is highest
cells_to_assign[[cluster_here]] <- intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% to_assign_cluster))[which(apply(mean_exprs[[cluster_here]], 1, function(x){if (x[1] >= max(x)) return(1) else {return(0)}}) == 1)]
pd_norm$cell_types_cl_all[cells_to_assign[[cluster_here]]] <- "epithelial"
table(pd_norm$cell_types_cl_all[which(pd_norm$Cluster == cluster_here)])
# cluster 8: try to assign the unknown and undecided cells to epithelial
cluster_here <- 8
table(pd_norm$cell_types_markers[which(pd_norm$Cluster == cluster_here)])
to_assign_cluster <- c("undecided", "unknown")
mean_exprs[[cluster_here]] <- compute_mean_expr_types(types = c("epithelial", "immune", "other"), mat_expr = mat_norm, 
                                                      cells_pos = intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% to_assign_cluster)),
                                                      epithelial_markers = epithelial_markers, immune_markers = immune_markers, other_markers = other_markers)
# only assign the cells whose epithelial expression is highest
cells_to_assign[[cluster_here]] <- intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% to_assign_cluster))[which(apply(mean_exprs[[cluster_here]], 1, function(x){if (x[1] >= max(x)) return(1) else {return(0)}}) == 1)]
pd_norm$cell_types_cl_all[cells_to_assign[[cluster_here]]] <- "epithelial"
table(pd_norm$cell_types_cl_all[which(pd_norm$Cluster == cluster_here)])
# also check whether the stroma cells have higher epithelial expression
mean_reassign_exprs[[cluster_here]] <- compute_mean_expr_types(types = c("epithelial", "immune", "other"), mat_expr = mat_norm, 
                                                               cells_pos = intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% "stroma")),
                                                               epithelial_markers = epithelial_markers, immune_markers = immune_markers, other_markers = other_markers)
# only re-assign the stroma cells if their epithelial expression is higher
cells_to_reassign[[cluster_here]] <- intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% "stroma"))[which(apply(mean_reassign_exprs[[cluster_here]], 1, function(x){if (x[1] >= max(x)) return(1) else {return(0)}}) == 1)]
pd_norm$cell_types_cl_all[cells_to_reassign[[cluster_here]]] <- "epithelial"
table(pd_norm$cell_types_cl_all[which(pd_norm$Cluster == cluster_here)])
# also check whether the Bcells or the Tcell have higher epithelial expression
mean_reassign_exprs[[cluster_here]] <- compute_mean_expr_types(types = c("epithelial", "immune", "other"), mat_expr = mat_norm, 
                                                               cells_pos = intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% c("Bcell", "Tcell"))),
                                                               epithelial_markers = epithelial_markers, immune_markers = immune_markers, other_markers = other_markers)
# only re-assign the Bcells or Tcell if their epithelial expression is higher
cells_to_reassign[[cluster_here]] <- intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% c("Bcell", "Tcell")))[which(apply(mean_reassign_exprs[[cluster_here]], 1, function(x){if (x[1] >= max(x)) return(1) else {return(0)}}) == 1)]
pd_norm$cell_types_cl_all[cells_to_reassign[[cluster_here]]] <- "epithelial"
table(pd_norm$cell_types_cl_all[which(pd_norm$Cluster == cluster_here)])
# cluster 9: try to assign the unknown and undecided cells to macrophage
cluster_here <- 9
table(pd_norm$cell_types_markers[which(pd_norm$Cluster == cluster_here)])
to_assign_cluster <- c("undecided", "unknown")
mean_exprs[[cluster_here]] <- compute_mean_expr_types(types = c("epithelial", "immune", "other"), mat_expr = mat_norm, 
                                                      cells_pos = intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% to_assign_cluster)),
                                                      epithelial_markers = epithelial_markers, immune_markers = immune_markers, other_markers = other_markers)
# only assign the cells whose mean immune expression is highest
cells_to_assign[[cluster_here]] <- intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% to_assign_cluster))[which(apply(mean_exprs[[cluster_here]], 1, function(x){if (x[2] >= max(x)) return(1) else {return(0)}}) == 1)]
pd_norm$cell_types_cl_all[cells_to_assign[[cluster_here]]] <- "macrophage"
table(pd_norm$cell_types_cl_all[which(pd_norm$Cluster == cluster_here)])
# also check whether the epithelial cells have high immune expression
mean_reassign_exprs[[cluster_here]] <- compute_mean_expr_types(types = c("epithelial", "immune", "other"), mat_expr = mat_norm, 
                                                               cells_pos = intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% "epithelial")),
                                                               epithelial_markers = epithelial_markers, immune_markers = immune_markers, other_markers = other_markers)
# only re-assign the epithelial cells whose mean immune expression is highest
cells_to_reassign[[cluster_here]] <- intersect(which(pd_norm$Cluster == cluster_here), which(pd_norm$cell_types_markers %in% "epithelial"))[which(apply(mean_reassign_exprs[[cluster_here]], 1, function(x){if (x[2] >= max(x)) return(1) else {return(0)}}) == 1)]
pd_norm$cell_types_cl_all[cells_to_reassign[[cluster_here]]] <- "macrophage"
table(pd_norm$cell_types_cl_all[which(pd_norm$Cluster == cluster_here)])

colData(sceset_final)$cell_types_cl_all <- pd_norm$cell_types_cl_all

## remove the unknown and undecided cells
unkund <- which(pd_norm$cell_types_cl_all %in% c("undecided", "unknown"))
# update all the matrices
sceset_ct <- sceset_final[,-unkund]
pd_ct <- colData(sceset_ct)
mat_ct <- assays(sceset_ct)$exprs
mats_ct <- list()
pds_ct <- list()
for (i in 1:length(patients_now)) {
  mats_ct[[i]] <- mat_ct[,pd_ct$patient == patients_now[i]]
  pds_ct[[i]] <- pd_ct[pd_ct$patient == patients_now[i],]
}
names(mats_ct) <- patients_now
names(pds_ct) <- patients_now
# output 
save(sceset_ct, file = "sceset_ct.RData")

## barplot cell types
match_celltype_levels <- c("epithelial", "stroma", "endothelial", "Tcell", "Bcell", "macrophage")
tbl_pd_ct <- as_tibble(pd_ct)
tbl_pd_ct <- tbl_pd_ct %>%
  group_by(patient) %>%
  mutate(cell_types_cl_all = factor(cell_types_cl_all, levels = match_celltype_levels)) %>%
  arrange(cell_types_cl_all)
# Figure 1C
ggplot() +
  geom_bar(data = tbl_pd_ct, aes(x = patient, fill = factor(cell_types_cl_all)), position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = anno_colors$tsne) +
  labs(fill = "cell type", y = "fraction of cells")

## barplot cell cycle phase for all patients
tbl_pd_cycle <- tbl_pd_ct %>%
  group_by(patient) %>%
  mutate(cycling_mel = factor(cycling_mel, levels = c("cycling", "non-cycling"))) %>%
  arrange(cycling_mel)
# Figure 1d
ggplot() +
  geom_bar(data = tbl_pd_cycle, aes(x = patient, fill = factor(cycling_mel)), position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = anno_colors$cycling) +
  labs(fill = "cycling status", y = "fraction of cells")
# end of Figure 1d

## patient cycling plots
p <- list()
for (i in 1:length(patients_now)) {
  #pdf(here("plots", paste("cycling_patient", patients_now[i], ".pdf", sep = "")))
  percent_epith <- length(intersect_all(which(pd_ct$patient == patients_now[i]), 
                                        which(pd_ct$cell_types_cl_all == "epithelial"), 
                                        which(pd_ct$cycling_mel == "cycling")))/length(intersect_all(
                                          which(pd_ct$patient == patients_now[i]), 
                                          which(pd_ct$cell_types_cl_all == "epithelial")))*100
  percent_all <- length(intersect_all(which(pd_ct$patient == patients_now[i]), 
                                      which(pd_ct$cycling_mel == "cycling")))/length(which(pd_ct$patient == patients_now[i]))*100
  if(i == 6)
    p[[i]] <- ggplot(as.data.frame(pds_ct[[i]]), aes(x = mel_scores_g1s, y = mel_scores_g2m)) +
    geom_rect(aes(xmin = median(pd_ct$mel_scores_g1s) + 2 * mad(pd_ct$mel_scores_g1s),
                  xmax = Inf, ymin = -Inf, ymax = Inf),
              fill = "gainsboro", alpha = 0.05) +
    geom_rect(aes(ymin = median(pd_ct$mel_scores_g2m) + 2 * mad(pd_ct$mel_scores_g2m),
                  ymax = Inf, xmin = -Inf, xmax = Inf),
              fill = "gainsboro", alpha = 0.05) +
    geom_point(aes(col = factor(cell_types_cl_all, levels = names(anno_colors$tsne)), 
                   shape = factor(cycling_mel)), size  = 5) +
    xlim(-0.15, 2) + ylim(-0.15, 2.8) +
    labs(col = "cell type", shape = "cycling", x = "G1S score", y = "G2M score", 
         title = paste("patient ", patients_now[i], " (", round(percent_all), "% cycling cells)", sep = "")) + 
    scale_color_manual(values = anno_colors$tsne) + theme(legend.position = "none", panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                          panel.background = element_blank(), axis.line = element_line(colour = "black"))
  
  if(i != 6)
    p[[i]] <- ggplot(as.data.frame(pds_ct[[i]]), aes(x = mel_scores_g1s, y = mel_scores_g2m)) +
    geom_rect(aes(xmin = median(pd_ct$mel_scores_g1s) + 2 * mad(pd_ct$mel_scores_g1s),
                  xmax = Inf, ymin = -Inf, ymax = Inf),
              fill = "gainsboro", alpha = 0.05) +
    geom_rect(aes(ymin = median(pd_ct$mel_scores_g2m) + 2 * mad(pd_ct$mel_scores_g2m),
                  ymax = Inf, xmin = -Inf, xmax = Inf),
              fill = "gainsboro", alpha = 0.05) +
    geom_point(aes(col = factor(cell_types_cl_all, levels = names(anno_colors$tsne)), 
                   shape = factor(cycling_mel)), size  = 5) +
    xlim(-0.15, 2) + ylim(-0.15, 2.8) +
    labs(col = "cell type", shape = "cycling", x = "G1S score", y = "G2M score", 
         title = paste("patient ", patients_now[i], " (", round(percent_all), "% cycling cells)", sep = "")) + 
    scale_color_manual(values = anno_colors$tsne) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                                          panel.background = element_blank(), axis.line = element_line(colour = "black"))
}
# Figure 1e
p[[6]] + p[[2]]
# end of Figure 1e

