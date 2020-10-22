# Yige Wu @WashU Sep 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 7
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input summary for cell-cell interactin
summary_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/summarize_interactions/rank_across_celltypes_by_pair_by_avgsigmean/20201012.v1/cellphonedb.summary_across_celltypes_by_pair.min5.tsv")
## input filtered summary data
interactions_filtered_df <- fread(data.table = F, input = "./Resources/snRNA_Processed_Data/CellPhoneDB/Filtered_Interactions/Values_averaged_mean_expr_of_all_interactions_stroma_immune_20200818_seen_7+_grouped_by_Inter.group_strongest_in_given_pair.txt")

# specify pairs to filter -------------------------------------------------
summary_df <- summary_df %>%
  mutate(paired_celltypes_group = paired_cellgroups.general)
summary_filtered_df <- summary_df %>%
  filter(pair_cell.types %in% interactions_filtered_df$pair_cell.types) %>%
  arrange(desc(paired_celltypes_group), desc(avg_sig_mean))
interactions_ordered_df <- interactions_filtered_df %>%
  filter(pair_cell.types %in% summary_filtered_df$pair_cell.types)
interacting_pairs_process <- rev(interactions_ordered_df$interacting_pair[!duplicated(interactions_ordered_df$interacting_pair)])
interacting_pairs_process <- summary_filtered_df$interacting_pair[!duplicated(summary_filtered_df$interacting_pair)]
interacting_pairs_process
length(interacting_pairs_process)
# paired_celltypes_process <- summary_filtered_df$paired_celltypes
paired_celltypes_process <- c("Endothelial cells|Macrophages", "Macrophages|Endothelial cells", "Macrophages|Myofibroblasts",
                              "Endothelial cells|Endothelial cells", "Endothelial cells|Myofibroblasts", "Fibroblasts|Endothelial cells", "Myofibroblasts|Endothelial cells", "Myofibroblasts|Fibroblasts", "Myofibroblasts|Myofibroblasts",
                              "Macrophages|Tumor cells", "TRM|Tumor cells", "Tumor cells|Macrophages", "Tumor cells|TRM",
                              "Endothelial cells|Tumor cells", "Myofibroblasts|Tumor cells", "Tumor cells|Endothelial cells", "Tumor cells|Fibroblasts", "Tumor cells|Myofibroblasts")
celltypes.source2target_process <- c("Endothelial cells->Tumor cells", "Myofibroblasts->Tumor cells", "Tumor cells->Endothelial cells", "Fibroblasts->Tumor cells", "Tumor cells->Myofibroblasts")
# make data matrix for the heatmap body -------------------------------------------
plot_data_df <- summary_df %>%
  filter(interacting_pair %in% interacting_pairs_process) %>%
  filter(paired_celltypes %in% paired_celltypes_process)
plot_data_wide_df <- dcast(data = plot_data_df, formula = interacting_pair ~ celltypes.source2target, value.var = "avg_sig_mean")
plot_data_mat <- as.matrix(plot_data_wide_df[, -1])
rownames(plot_data_mat) <- plot_data_wide_df$interacting_pair
plot_data_mat <- plot_data_mat[interacting_pairs_process,]
plot_data_mat[is.na(plot_data_mat)] <- -1

# get ids -----------------------------------------------------------------
interacting_pair_plot <- rownames(plot_data_mat)
celltypes.source2target_plot <- colnames(plot_data_mat)

# specify colors ----------------------------------------------------------
color_red <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
color_orange <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[5]
color_yellow <- RColorBrewer::brewer.pal(n = 9, name = "YlGnBu")[1]
## specify color for NA values
color_na <- "grey70"
## make color function for heatmap body colors
summary(as.vector(plot_data_mat[plot_data_mat != -1]))
colors_heatmapbody <- colorRamp2(c(-1, 0, 2, 4, 6, 8, 10, 12, 14, 16), 
                                 c(color_na, brewer.pal(n = 9, name = "YlOrRd")))
colors_heatmapbody_legand <- colorRamp2(c(0, 2, 4, 6, 8, 10, 12, 14, 16), 
                                        c(brewer.pal(n = 9, name = "YlOrRd")))
colors_heatmapbody <- colorRamp2(c(-1, 0, 2, 4, 6, 8, 10, 12), 
                                 c(color_na, brewer.pal(n = 7, name = "YlOrRd")))
colors_heatmapbody_legand <- colorRamp2(c(0, 2, 4, 6, 8, 10, 12), 
                                        c(brewer.pal(n = 7, name = "YlOrRd")))
## cell type colors
colors_celltype <- RColorBrewer::brewer.pal(n = 6, name = "Set1")
names(colors_celltype) <- c("Tumor cells", "Endothelial cells", "Myofibroblasts", "Fibroblasts", "Macrophages", "TRM")
## make colors for pathways
colors_pathway <- RColorBrewer::brewer.pal(n = 10, name = "Set3")
names(colors_pathway) <- c("VEGFR", "EGFR", "ERBB", "Ephrin", "FGFR", "WNT", "EMT", "TGFb", "", "Integrin")

# make row annotation -----------------------------------------------------
## the number of samples with significant strength
number_sig_samples_vec <- mapvalues(x = interacting_pair_plot, from = summary_filtered_df$interacting_pair, to = as.vector(summary_filtered_df$number_sig_cases))
number_sig_samples_vec <- as.numeric(number_sig_samples_vec)
number_sig_samples_vec
## make text
gene_source_vec <- mapvalues(x = interacting_pair_plot, from = summary_filtered_df$interacting_pair, to = as.vector(summary_filtered_df$gene.source))
gene_target_vec <- mapvalues(x = interacting_pair_plot, from = summary_filtered_df$interacting_pair, to = as.vector(summary_filtered_df$gene.target))
## make cell types
celltype_source_vec <- mapvalues(x = interacting_pair_plot, from = summary_filtered_df$interacting_pair, to = as.vector(summary_filtered_df$Cell_type.source))
celltype_target_vec <- mapvalues(x = interacting_pair_plot, from = summary_filtered_df$interacting_pair, to = as.vector(summary_filtered_df$Cell_type.target))
## make interaction direciton
# make annotation object
rowanno_obj1 <- rowAnnotation(
  # Celltype_ligand = anno_simple(x = vector(mode = "numeric", length = length(celltype_source_vec)), 
  #                               gp = gpar(fill = "white"),
  #                               pch = 21, pt_gp = gpar(fill = colors_celltype[celltype_source_vec])),
  Gene_source = anno_text(x = gene_source_vec, 
                          gp = gpar(fontface = "italic", fontsize = 30, border = "black"), 
                          location = 0.5,
                          just = "center"),
  # Celltype_receptor = anno_simple(x = vector(mode = "numeric", length = length(celltype_target_vec)), 
  #                                 gp = gpar(fill = "white"),
  #                                 pch = 21, pt_gp = gpar(fill = colors_celltype[celltype_target_vec])),
  Gene_target = anno_text(x = gene_target_vec, 
                          gp = gpar(fontface = "italic", fontsize = 30, border = "black"), 
                          location = 0.5,
                          just = "center"),
  # Pathway = anno_text(x = pathway_vec, gp = gpar(fill = colors_pathway[pathway_vec])),
  annotation_name_side = "top")
rowanno_obj2 <- rowAnnotation(Number_samples = anno_barplot(x = number_sig_samples_vec, width = unit(x = 2, "cm")),
                              annotation_name_side = "bottom")


# make column annotation -----------------------------------------------------
paired_cellgroups_detailed <- mapvalues(x = celltypes.source2target_plot, from = summary_df$celltypes.source2target, to = as.vector(summary_df$paired_celltypes_group)) 
col_celltypes.source_vec <- mapvalues(x = celltypes.source2target_plot, from = summary_df$celltypes.source2target, to = as.vector(summary_df$Cell_type.source))
col_celltypes.target_vec <- mapvalues(x = celltypes.source2target_plot, from = summary_df$celltypes.source2target, to = as.vector(summary_df$Cell_type.target))
## get direction
direction_pch_vec <- 25
colanno_obj <- HeatmapAnnotation(Celltype.Ligand = anno_simple(x = vector(mode = "numeric", length = length(col_celltypes.source_vec)), 
                                                               gp = gpar(fill = "white"), height = unit(10, "mm"), 
                                                               border = F,
                                                               pch = 21, pt_size = unit(10, "mm"), pt_gp = gpar(fill = colors_celltype[col_celltypes.source_vec])), 
                                 Direction = anno_simple(x = vector(mode = "numeric", length = length(col_celltypes.source_vec)), 
                                                         gp = gpar(fill = "white"), height = unit(8, "mm"), 
                                                         pch = direction_pch_vec, pt_size = unit(8, "mm"), pt_gp = gpar(fill = "black")),
                                 Celltype.Receptor = anno_simple(x = vector(mode = "numeric", length = length(col_celltypes.source_vec)), 
                                                                 gp = gpar(fill = "white"), height = unit(10, "mm"), 
                                                                 pch = 21, pt_size = unit(10, "mm"), pt_gp = gpar(fill = colors_celltype[col_celltypes.target_vec])),
                                 annotation_name_side = "right", annotation_name_gp = gpar(fontsize = 25))


# make column split -------------------------------------------------------
column_split_vec <- paired_cellgroups_detailed
column_split_factor <- factor(x = column_split_vec, levels = c("Stroma&Immune", "Stroma&Stroma", "Tumor&Immune", "Tumor&Stroma"))

# plot hetamp body --------------------------------------------------------
plot_size_df <- dcast(data = plot_data_df, formula = interacting_pair ~ celltypes.source2target, value.var = "number_sig_cases")
plot_size_mat <- as.matrix(plot_size_df[, -1])
rownames(plot_size_mat) <- plot_size_df$interacting_pair
plot_size_mat <- plot_size_mat[interacting_pairs_process,]

p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat,
                             col = colors_heatmapbody,
                             width = unit(ncol(plot_data_mat), "cm"), height = unit(nrow(plot_data_mat), "cm"),
                             border = "black",
                             ## cell function
                             rect_gp = gpar(type = "none"), 
                             cell_fun = function(j, i, x, y, width, height, fill) {
                               grid.rect(x = x, y = y, width = width, height = height, 
                                         gp = gpar(col = color_na, fill = NA))
                               if (plot_data_mat[i, j] > 0) {
                                 # grid.circle(x = x, y = y, r = ((plot_size_mat[i, j])^(1/3)+0.5) * min(unit.c(width, height)), 
                                 #             gp = gpar(fill = colors_heatmapbody(plot_data_mat[i, j]), col = "black"))
                                 grid.rect(x = x, y = y, 
                                           width = (sqrt(plot_size_mat[i, j])/sqrt(30)) * min(unit.c(width)),
                                           height = (sqrt(plot_size_mat[i, j])/sqrt(30)) * min(unit.c(height)),
                                             # width = ((plot_size_mat[i, j]^(1/3))/2) * min(unit.c(width)),
                                             # height = ((plot_size_mat[i, j]^(1/3))/2) * min(unit.c(height)),
                                             gp = gpar(fill = colors_heatmapbody(plot_data_mat[i, j]), col = "black"))
                               }},
                             ## row parameters
                             
                             # row_names_side = "left",
                             # row_names_max_width = unit(10, "cm"),
                             show_row_dend = F, cluster_rows = F,
                             left_annotation = rowanno_obj1, 
                             # right_annotation = rowanno_obj2,
                             show_row_names = F,
                             
                             ## column parameters
                             column_gap = unit(x = 0, units = "cm"),
                             column_split = column_split_factor, cluster_column_slices = F, cluster_columns = F,
                             top_annotation = colanno_obj, show_column_names = F,
                             column_names_side = "top", column_title_rot = 90, column_title_gp = gpar(fontsize = 30),
                             show_column_dend = F,
                             show_heatmap_legend = F)
p

# make legend list --------------------------------------------------------
list_lgd = list(
  Legend(col_fun = colors_heatmapbody_legand, 
         title = "Interaction score", 
         labels_gp = gpar(fontsize=25),
         title_gp = gpar(fontsize = 30),
         legend_width = unit(5, "cm"),
         legend_height = unit(5, "cm"),
         direction = "horizontal"),
  Legend(title = NULL, 
         # title_gp = gpar(fontsize = 15),
         labels = "Not significant",
         labels_gp = gpar(fontsize=25), 
         legend_gp = gpar(fill = color_na),
         grid_height = unit(1, "cm"),
         direction = "vertical"),
  Legend(labels = names(colors_celltype),
         title = "Cell type",
         legend_gp = gpar(fill = colors_celltype), 
         labels_gp = gpar(fontsize=25), title_gp = gpar(fontsize = 30),
         grid_height = unit(1, "cm"),
         grid_width = unit(1, "cm")))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "interactions", ".png")
png(file2write, width = 2600, height = 4000, res = 150)
draw(object = p, 
     annotation_legend_side = "left", annotation_legend_list = list_lgd)
dev.off()

file2write <- paste0(dir_out, "interactions", ".pdf")
pdf(file2write, width = 18, height = 25)
draw(object = p, 
     annotation_legend_side = "left", annotation_legend_list = list_lgd)
dev.off()

