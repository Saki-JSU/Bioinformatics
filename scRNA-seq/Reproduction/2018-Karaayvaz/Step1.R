# Figure 1b
# initial environment
rm(list=ls())
gc()
setwd("C:/Users/saki/Desktop/2018-Karaayvaz/code")

# library
library(SingleCellExperiment)
library(ComplexHeatmap)
library(circlize)
library(broman)

# load functions
source("funcs.R")
source("funcs_markers.R")

## read in normalized data and phenotypic information
mat_norm <- read.table("../data/norm_data.txt", sep = "\t", header = TRUE)
pd_norm <- read.table("../data/pd_norm.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
fd_norm <- read.table("../data/fd_norm.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
sceset_final <- SingleCellExperiment(assays = list(exprs = as.matrix(mat_norm)),
                                     colData = pd_norm, rowData = fd_norm)

## colors for plotting
epithelial_col <- brocolors("crayons")["Maroon"]
basal_epithelial_col <- brocolors("crayons")["Red"]
luminal_epithelial_col <- brocolors("crayons")["Sunset Orange"]
luminal_progenitor_col <- brocolors("crayons")["Salmon"]

stroma_col <- brocolors("crayons")["Aquamarine"]
endothelial_col <- brocolors("crayons")["Wisteria"]

PTPRC_col <- brocolors("crayons")["Inchworm"]
t_cell_col <- brocolors("crayons")["Screamin' Green"]
b_cell_col <- brocolors("crayons")["Fern"]
macrophage_col <- brocolors("crayons")["Tropical Rain Forest"]

marker_cols <- c("epithelial" = unname(epithelial_col), "basal epithelial" = unname(basal_epithelial_col), 
                 "luminal epithelial" = unname(luminal_epithelial_col), "luminal progenitor" = unname(luminal_progenitor_col),
                 "stroma" = unname(stroma_col), "endothelial" = unname(endothelial_col), 
                 "immune" = unname(PTPRC_col), "T cell" = unname(t_cell_col), "B cell" = unname(b_cell_col), 
                 "macrophage" = unname(macrophage_col))
cycling_mel_cols <- c("non-cycling" = "gainsboro",
                      "cycling" = unname(brocolors("crayons")["Mulberry"]))
depletion_cols <- c("depleted" = unname(brocolors("crayons")["White"]), "not depleted" = unname(brocolors("crayons")["Red"]))
pats_cols <- c("PT039" = unname(brocolors("crayons")["Orange Red"]), "PT058" = unname(brocolors("crayons")["Orange"]), 
               "PT081" = unname(brocolors("crayons")["Pink Flamingo"]), "PT084" = unname(brocolors("crayons")["Fern"]), 
               "PT089" = unname(brocolors("crayons")["Blue Violet"]), "PT126" = unname(brocolors("crayons")["Sky Blue"]))
tsne_cols <- c("epithelial" = unname(basal_epithelial_col), "stroma" = unname(stroma_col), "endothelial" = unname(endothelial_col),
               "Tcell" = unname(t_cell_col), "Bcell" = unname(b_cell_col), "macrophage" = unname(macrophage_col))
#tsne_cols_unknown <- c(tsne_cols, "unknown" = "black", "undecided" = "black")
anno_colors <- list("marker" = marker_cols, "cycling" = cycling_mel_cols, "immune depletion" = depletion_cols, 
                    "patient" = pats_cols, "tsne" = tsne_cols)
# output
save(anno_colors, file = "anno_colors.RData")

## cycling cells identification
## cell cycle assignment as done in melanoma in Tirosch et al 2016
# g1s and g2m periods
melanoma_cellcycle <- read.table("../data/melanoma_cellcycle.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
melanoma_g1s <- melanoma_cellcycle$G1S
melanoma_g1s <- match_clean_vector_genes(mat_norm, melanoma_g1s)
scores_g1s <- avg_expr_genes(mat_norm, melanoma_g1s$index)
melanoma_g2m <- melanoma_cellcycle$G2M
melanoma_g2m <- match_clean_vector_genes(mat_norm, melanoma_g2m)
scores_g2m <- avg_expr_genes(mat_norm, melanoma_g2m$index)
# classification
cycling_mel <- rep(NA, length(scores_g1s))
cycling_mel[scores_g1s >= (median(scores_g1s) + 2 * mad(scores_g1s)) & scores_g2m < (median(scores_g2m) + 2 * mad(scores_g2m))] <- "cycling"
cycling_mel[scores_g1s < (median(scores_g1s) + 2 * mad(scores_g1s)) & scores_g2m >= (median(scores_g2m) + 2 * mad(scores_g2m))] <- "cycling"
cycling_mel[scores_g1s < (median(scores_g1s) + 2 * mad(scores_g1s)) & scores_g2m < (median(scores_g2m) + 2 * mad(scores_g2m))] <- "non-cycling"
cycling_mel[scores_g1s >= (median(scores_g1s) + 2 * mad(scores_g1s)) & scores_g2m >= (median(scores_g2m) + 2 * mad(scores_g2m))] <- "cycling"
# update colData and pd_norm
colData(sceset_final)$mel_scores_g1s <- scores_g1s
colData(sceset_final)$mel_scores_g2m <- scores_g2m
colData(sceset_final)$cycling_mel <- cycling_mel
pd_norm <- colData(sceset_final)
pd_norm <- as.data.frame(pd_norm)

## individual patient
patients_now <- c()
mats_now <- list()
pds_now <- list()
for (i in 1:length(unique(pd_norm$patient))) {
  patients_now[i] <- sort(unique(pd_norm$patient))[i]
  mats_now[[i]] <- mat_norm[, pd_norm$patient == patients_now[i]]
  pds_now[[i]] <- pd_norm[pd_norm$patient == patients_now[i],]
}
names(mats_now) <- patients_now
names(pds_now) <- patients_now
# output
save(sceset_final, file = "sceset_final.RData")

## cell types identification (via markers and clustering)
all_markers <- read.table("../data/markers_clean.txt", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
markers <- unique(all_markers[all_markers$gene %in% rowData(sceset_final)$hgnc_symbol, ])
# output
save(markers, file = "markers.RData")

# marker names
# small type
markers$type_heatmap <- markers$type
markers$type_heatmap[which(markers$type_heatmap == "luminalprogenitor")] <- "luminal progenitor"
markers$type_heatmap[which(markers$type_heatmap == "luminalepithelial")] <- "luminal epithelial"
markers$type_heatmap[which(markers$type_heatmap == "basalepithelial")] <- "basal epithelial"
markers$type_heatmap[which(markers$type_heatmap %in% c("EPCAM", "EGFR", "CDH1"))] <- "epithelial"
markers$type_heatmap[which(markers$type_heatmap == "Bcell")] <- "B cell"
markers$type_heatmap[which(markers$type_heatmap == "Tcell")] <- "T cell"
# large type
markers$type_long_heatmap <- markers$type_long
markers$type_long_heatmap[which(markers$type == "stroma")] <- "stroma"
markers$type_long_heatmap[which(markers$type == "endothelial")] <- "endothelial"

# Figure 1b
# colors for markers name text on the right
colors_markers_ch <- markers$type_heatmap
for (i in c(1:length(names(anno_colors$marker)))) {
  colors_markers_ch <- replace(colors_markers_ch, colors_markers_ch == names(anno_colors$marker)[i], anno_colors$marker[i])
}
# change the order of large type text on the left
splits_ch <- as.factor(markers$type_long_heatmap)
splits_ch <- factor(splits_ch, levels(splits_ch)[c(2,4,1,3)])
# annotation on the left
colors_anno_markers_ch <- as.factor(markers$type_heatmap)
colors_anno_markers_ch <- factor(colors_anno_markers_ch, levels(colors_anno_markers_ch)[c(4,2,5,6,9,3,8,10,1,7)])
ha_rows <- HeatmapAnnotation(df = data.frame(type = colors_anno_markers_ch),
                             annotation_legend_param = list(type = list(ncol = 2, title = "cell type", title_position = "topcenter")),
                             which = "row", col = list("type" = anno_colors$marker), annotation_width = unit(3, "mm"),
                             show_annotation_name = FALSE)
# annotation on the up
cycling_now <- list()
ha_cols_up <- list()
for(i in 1:length(patients_now)){
  # cycling status
  cycling_now[[i]] <- pds_now[[i]][,"cycling_mel"]
  # first subfig with annotation
  if (i == 1)
    ha_cols_up[[i]] <- HeatmapAnnotation(df = data.frame(cycling = cycling_now[[i]]), 
                                         col = list(cycling = anno_colors$cycling), 
                                         show_annotation_name = TRUE, 
                                         annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 11),
                                         annotation_legend_param = list(list(title_position = "topcenter",
                                                                             title = c("cycling status"))),
                                         gap = unit(c(1, 1), "mm"))
  # other subfig without annotations
  if (i > 1)
    ha_cols_up[[i]] <- HeatmapAnnotation(df = data.frame(cycling = cycling_now[[i]]), 
                                         col = list(cycling = anno_colors$cycling), 
                                         show_legend = FALSE,
                                         gap = unit(c(1, 1), "mm"))
}

# annotation on the bottom
depletion_now <- list()
ha_cols_bottom <- list()
for(i in 1:length(patients_now)){
  # cycling status
  depletion_now[[i]] <- pds_now[[i]][,"depletion_batch"]
  # first subfig with annotation
  if (i == 1)
    ha_cols_bottom[[i]] <- HeatmapAnnotation(df = data.frame('CD45' = depletion_now[[i]]), 
                                             col = list('CD45' = c("depleted_yes" = "gainsboro", "depleted_no" = "gray54")),
                                             show_annotation_name = TRUE,
                                             annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 11),
                                             annotation_legend_param = list(title = "CD45 status", 
                                                                            title_position = "topcenter", 
                                                                            at = c("depleted_yes", "depleted_no"),
                                                                            labels = c("CD45 depleted","CD45 unselected")),
                                             show_legend = TRUE,
                                             gap = unit(c(1), "mm"))
  
  if (i > 1)
    ha_cols_bottom[[i]] <- HeatmapAnnotation(df = data.frame('CD45' = depletion_now[[i]]), 
                                             col = list('CD45' = c("depleted_yes" = "gainsboro", "depleted_no" = "gray54")),
                                             show_annotation_name = FALSE,
                                             show_legend = FALSE,
                                             gap = unit(c(1), "mm"))
}
names(cycling_now) <- patients_now
names(depletion_now) <- patients_now
names(ha_cols_up) <- patients_now
names(ha_cols_bottom) <- patients_now
# plot figures
ht_list <- ha_rows + 
  Heatmap(mats_now[[1]][match(markers$gene,rownames(mats_now[[1]])),], 
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE,
          name = patients_now[1], clustering_distance_columns = "euclidean", row_names_side = "left", 
          row_names_gp = gpar(fontsize = 10, col = colors_markers_ch),
          split = splits_ch, gap = unit(1, "mm"), column_title = patients_now[1], 
          column_title_gp = gpar(fontsize = 11),
          row_title_gp = gpar(font = 11), top_annotation = ha_cols_up[[1]], 
          heatmap_legend_param = list(title = "expression", title_position = "topcenter", color_bar = "continuous"),
          bottom_annotation = ha_cols_bottom[[1]],
          col = colorRamp2(c(-1, 3.5, 8), c("blue", "white", "red"))) + 
  Heatmap(mats_now[[2]][match(markers$gene,rownames(mats_now[[1]])),],
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE,
          name = patients_now[2], clustering_distance_columns = "euclidean",
          show_row_names = FALSE, column_title = patients_now[2],
          column_title_gp = gpar(fontsize = 11),
          show_heatmap_legend = FALSE,
          row_title_gp = gpar(font = 11), top_annotation = ha_cols_up[[2]], bottom_annotation = ha_cols_bottom[[2]],
          gap = unit(1, "mm"),
          col = colorRamp2(c(-1, 3.5, 8), c("blue", "white", "red"))) +
  Heatmap(mats_now[[3]][match(markers$gene,rownames(mats_now[[3]])),],
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE,
          name = patients_now[3], clustering_distance_columns = "euclidean",
          show_row_names = FALSE, column_title = patients_now[3],
          column_title_gp = gpar(fontsize = 11),
          show_heatmap_legend = FALSE,
          row_title_gp = gpar(font = 11), top_annotation = ha_cols_up[[3]], bottom_annotation = ha_cols_bottom[[3]],
          gap = unit(1, "mm"),
          col = colorRamp2(c(-1, 3.5, 8), c("blue", "white", "red"))) +
  Heatmap(mats_now[[4]][match(markers$gene,rownames(mats_now[[4]])),],
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE,
          name = patients_now[4], clustering_distance_columns = "euclidean",
          show_row_names = FALSE, column_title = patients_now[4], 
          column_title_gp = gpar(fontsize = 11),
          show_heatmap_legend = FALSE,
          row_title_gp = gpar(font = 11), top_annotation = ha_cols_up[[4]], bottom_annotation = ha_cols_bottom[[4]],
          gap = unit(1, "mm"),
          col = colorRamp2(c(-1, 3.5, 8), c("blue", "white", "red"))) +
  Heatmap(mats_now[[5]][match(markers$gene,rownames(mats_now[[5]])),], 
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE,
          name = patients_now[5], clustering_distance_columns = "euclidean",
          show_row_names = FALSE, column_title = patients_now[5], 
          column_title_gp = gpar(fontsize = 11), 
          show_heatmap_legend = FALSE,
          row_title_gp = gpar(font = 11), top_annotation = ha_cols_up[[5]], bottom_annotation = ha_cols_bottom[[5]],
          gap = unit(1, "mm"),
          col = colorRamp2(c(-1, 3.5, 8), c("blue", "white", "red"))) + 
  Heatmap(mats_now[[6]][match(markers$gene,rownames(mats_now[[6]])),], 
          cluster_rows = FALSE, show_column_dend = TRUE, show_column_names = FALSE,
          name = patients_now[6], clustering_distance_columns = "euclidean",
          row_names_gp = gpar(fontsize = 10, col = colors_markers_ch), 
          split = splits_ch, column_title = patients_now[6], 
          column_title_gp = gpar(fontsize = 11),
          row_title_gp = gpar(font = 11), top_annotation = ha_cols_up[[6]], bottom_annotation = ha_cols_bottom[[6]],
          show_heatmap_legend = FALSE,
          gap = unit(1, "mm"),col = colorRamp2(c(-1, 3.5, 8), c("blue", "white", "red")))
draw(ht_list, gap = unit(0.1, "cm"), heatmap_legend_side = "bottom", annotation_legend_side = "bottom")
# end of figure 1b
