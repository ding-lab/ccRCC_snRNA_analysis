# Yige Wu @WashU Sep 2020
## 2020-11-16 Ref: https://www.ncbi.nlm.nih.gov/books/NBK12423/

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input summary for cell-cell interactin
summary_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/summarize_interactions/rank_across_celltypes_by_pair_by_avgsigmean/20201012.v1/cellphonedb.summary_across_celltypes_by_pair.min5.tsv")
## input filtered summary data
# interactions_filtered_df1 <- fread(data.table = F, input = "./Resources/snRNA_Processed_Data/CellPhoneDB/Filtered_Interactions/Values_averaged_mean_expr_of_all_interactions_stroma_immune_20200818_seen_7+_grouped_by_Inter.group_strongest_in_given_pair.txt")
interactions_filtered_df1 <- fread(data.table = F, input = "./Data_Freezes/V1/snRNA/Cell_Cell_Interactions/Values_by_sample_mean_expr_of_curated_interactions_top3ct_stroma_immune_20200818_seen_5+_grouped_by_Inter.group_strongest_in_given_pair.txt")

# specify pairs to filter -------------------------------------------------
summary_df <- summary_df %>%
  mutate(paired_celltypes_group = paired_cellgroups.general)
interactions_filtered_df2 <- interactions_filtered_df1 %>%
  filter(Inter.group != "Immune:Stroma")

# make data matrix for the heatmap body -------------------------------------------
# plot_data_df <- summary_df %>%
#   filter(interacting_pair %in% interacting_pairs_process) %>%
#   filter(paired_celltypes %in% paired_celltypes_process)
plot_data_df <- summary_df %>%
  filter(pair_cell.types %in% interactions_filtered_df2$pair_cell.types) %>%
  filter(grepl(x = gene.target, pattern = "PDGFR|CSF1R|KIT|FLT3|FLT1|KDR|FLT4|EGFR|ERBB|FGFR|INSR|IGF1R|INSRR|MET|MST1R|NTRK|AXL|TYRO3|MERTK|EPHA|EPHB|MUSK|RET|TEK|DDR|ROR|ROS1|RYK")) %>%
  filter(gene.target != "FLT1 complex") %>%
  filter(!(Cell_type.target %in% c("CD4 T-cells", "CD4 CTL"))) %>%
  mutate(interacting_pair_directed = paste0(gene.source, "_", gene.target)) %>%
  arrange(desc(paired_celltypes_group), desc(avg_sig_mean))
file2write <- paste0(dir_out, "ccRCC_sn_interactions_selected.", run_id, ".tsv")
write.table(x = plot_data_df, file = file2write, sep = "\t", row.names = F, quote = F)
plot_data_wide_df <- dcast(data = plot_data_df, formula = celltypes.source2target ~ interacting_pair_directed, value.var = "avg_sig_mean")
plot_data_mat <- as.matrix(plot_data_wide_df[, -1])
rownames(plot_data_mat) <- plot_data_wide_df$celltypes.source2target
plot_data_mat[is.na(plot_data_mat)] <- -1

# get ids including sorted-----------------------------------------------------------------
interacting_pair_plot <- colnames(plot_data_mat)
interacting_pair_plot_sorted <- c(
  ## tumor&stroma
  "VEGFA_FLT1", "VEGFA_KDR", "VEGFB_FLT1", "PGF_FLT1", "EFNA5_EPHA2", "EFNA5_EPHA4", 
  "EFNB2_EPHB1", "EFNB2_EPHA4", "EFNA1_EPHA7", "FGF2_FGFR2", "FGF2_FGFR1", "EFNA1_EPHA4", "NRG2_ERBB3", "FGF2_FGFR3",
  "HGF_MET", "FGF1_FGFR1",
  "EFNA5_EPHA3", "NRG3_ERBB4", "EFNA1_EPHA3", "NRG1_ERBB4", 
  "TGFA_EGFR", "EFNA5_EPHA7", "NRG1_ERBB3",
  ## tumor&tumor
  "TGFA_EGFR", "EFNA5_EPHA7", "EGF_EGFR", "EFNA5_EPHA4", "EFNA1_EPHA7", "NRG1_ERBB3", "FGF2_FGFR1", "NCAM1_FGFR1",
  ## stroma & stroma
  "VEGFA_FLT1", "PGF_FLT1", "EFNB2_EPHA4", "EFNB2_EPHB1", "EFNB2_EPHB4", "ANGPT2_TEK", "VEGFC_KDR", "VEGFC_FLT4", "EFNB1_EPHA4", "EFNA1_EPHA4", "FGF2_FGFR1", "EFNB1_EPHB4", "EFNA1_EPHA2",
  "VEGFB_FLT1", "ANGPT1_TEK", "FGF1_FGFR1",
  "EFNA5_EPHA7",
  "PDGFD_PDGFRB", "PDGFB_PDGFRB", "EFNA1_EPHA3", "EFNB2_EPHB6",
  "EFNA5_EPHA7", "TGFA_EGFR", 
  "EFNA5_EPHA3",
  "NRG3_ERBB4", "NTF3_NTRK3",
  ## tumor&immune
  "HBEGF_EGFR", "TGFA_EGFR", "IGF1_IGF1R", "HGF_MET",
  "IL34_CSF1R"
  )
interacting_pair_plot_sorted <- interacting_pair_plot_sorted[!duplicated(interacting_pair_plot_sorted)]
plot_data_mat <- plot_data_mat[,interacting_pair_plot_sorted]
celltypes.source2target_plot <- rownames(plot_data_mat)
celltypes.source2target_plot_sorted <- c("Tumor cells->Endothelial cells", "Endothelial cells->Tumor cells", 
                                         "Myofibroblasts->Tumor cells", "Tumor cells->Myofibroblasts", 
                                         "Tumor cells->Fibroblasts", "Fibroblasts->Tumor cells",
                                         "Tumor cells->Tumor cells",
                                         "Endothelial cells->Endothelial cells", "Myofibroblasts->Endothelial cells", 
                                         "Endothelial cells->Myofibroblasts", 
                                         "Myofibroblasts->Myofibroblasts", 
                                         "Fibroblasts->Endothelial cells", 
                                         "Fibroblasts->Fibroblasts", "Fibroblasts->Myofibroblasts",
                                         # "Endothelial cells->Fibroblasts",
                                         "Macrophages->Tumor cells", "Macrophages proliferating->Tumor cells", "TRM->Tumor cells",
                                         "Tumor cells->Macrophages", "Tumor cells->Macrophages proliferating",
                                         "NK cells weak->Tumor cells", "NK cells strong->Tumor cells"
)
plot_data_mat <- plot_data_mat[celltypes.source2target_plot_sorted, ]

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
unique(c(plot_data_df$Cell_type.source, plot_data_df$Cell_type.target))
colors_celltype <- rep(x = colors_cellgroup14[c("Tumor cells", "Endothelial cells", "Myofibroblasts", "Fibroblasts", "Macrophages", "NK cells")], c(1, 1, 1, 1, 3, 2))
names(colors_celltype) <- c("Tumor cells", "Endothelial cells", "Myofibroblasts", "Fibroblasts", 
                            "Macrophages", "Macrophages proliferating", "TRM", 
                            "NK cells strong", "NK cells weak")
# colors_celltype <- RColorBrewer::brewer.pal(n = 6, name = "Set1")
# names(colors_celltype) <- c("Tumor cells", "Endothelial cells", "Myofibroblasts", "Fibroblasts", "Macrophages", "TRM")
## make colors for pathways

# make column annotation -----------------------------------------------------
## make text
gene_source_vec <- mapvalues(x = interacting_pair_plot_sorted, from = plot_data_df$interacting_pair_directed, to = as.vector(plot_data_df$gene.source))
gene_target_vec <- mapvalues(x = interacting_pair_plot_sorted, from = plot_data_df$interacting_pair_directed, to = as.vector(plot_data_df$gene.target))
## make interaction direciton
# make annotation object
colanno_obj <- HeatmapAnnotation(
  Gene_source = anno_text(x = gene_source_vec,
                          gp = gpar(fontface = "italic", fontsize = 30, border = NA),
                          location = 0.5,
                          just = "center"),
  # Text_arrow = anno_text(x = rep(x = "->", length(gene_source_vec)),
  #                        gp = gpar(fontsize = 20, border = NA),
  #                        location = 0.5,
  #                        just = "center"),
  Gene_target = anno_text(x = gene_target_vec,
                          gp = gpar(fontface = "italic", fontsize = 30, border = NA),
                          location = 0.5,
                          just = "center"),
  annotation_name_side = "left")

# make row annotation -----------------------------------------------------
paired_cellgroups_detailed <- mapvalues(x = celltypes.source2target_plot_sorted, from = summary_df$celltypes.source2target, to = as.vector(summary_df$paired_celltypes_group)) 
col_celltypes.source_vec <- mapvalues(x = celltypes.source2target_plot_sorted, from = summary_df$celltypes.source2target, to = as.vector(summary_df$Cell_type.source))
col_celltypes.target_vec <- mapvalues(x = celltypes.source2target_plot_sorted, from = summary_df$celltypes.source2target, to = as.vector(summary_df$Cell_type.target))
## get direction
rowanno_obj <- rowAnnotation(Celltype.Ligand = anno_simple(x = vector(mode = "numeric", length = length(col_celltypes.source_vec)), 
                                                               gp = gpar(fill = NA, col = NA), width = unit(10, "mm"), 
                                                               pch = 21, 
                                                               pt_size = unit(10, "mm"), 
                                                               pt_gp = gpar(fill = colors_celltype[col_celltypes.source_vec], col = NA)), 
                                 # Text_arrow = anno_text(x = rep(x = "<-", length(col_celltypes.source_vec)),
                                 #                        gp = gpar(fontsize = 25, border = NA),
                                 #                        location = 0.5,
                                 #                        just = "center"),
                                 Celltype.Receptor = anno_simple(x = vector(mode = "numeric", length = length(col_celltypes.source_vec)),
                                                                 gp = gpar(fill = NA, col = NA), width = unit(10, "mm"), 
                                                                 pch = 21, 
                                                                 pt_size = unit(10, "mm"), 
                                                                 pt_gp = gpar(fill = colors_celltype[col_celltypes.target_vec], col = NA)),
                                 annotation_name_side = "bottom", annotation_name_gp = gpar(fontsize = 25))


# make column split -------------------------------------------------------
row_split_vec <- paired_cellgroups_detailed
row_split_factor <- factor(x = row_split_vec, levels = c("Tumor&Stroma", "Tumor&Tumor", "Stroma&Stroma", "Tumor&Immune", "Stroma&Immune"))

# plot hetamp body --------------------------------------------------------
plot_size_df <- dcast(data = plot_data_df, formula =celltypes.source2target ~ interacting_pair_directed, value.var = "number_sig_cases")
plot_size_mat <- as.matrix(plot_size_df[, -1])
rownames(plot_size_mat) <- plot_size_df[,1]
plot_size_mat <- plot_size_mat[celltypes.source2target_plot_sorted,interacting_pair_plot_sorted]

p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat,
                             col = colors_heatmapbody,
                             width = unit(ncol(plot_data_mat), "cm"), height = unit(nrow(plot_data_mat), "cm"),
                             # border = "black",
                             ## cell function
                             rect_gp = gpar(type = "none"), 
                             cell_fun = function(j, i, x, y, width, height, fill) {
                               grid.rect(x = x, y = y, width = width, height = height, 
                                         gp = gpar(col = "white", fill = "grey90"))
                               if (plot_data_mat[i, j] > 0) {
                                 # grid.circle(x = x, y = y, r = ((plot_size_mat[i, j])^(1/3)+0.5) * min(unit.c(width, height)), 
                                 #             gp = gpar(fill = colors_heatmapbody(plot_data_mat[i, j]), col = "black"))
                                 grid.rect(x = x, y = y, 
                                           width = (sqrt(plot_size_mat[i, j])/sqrt(30)) * min(unit.c(width)),
                                           height = (sqrt(plot_size_mat[i, j])/sqrt(30)) * min(unit.c(height)),
                                           # width = ((plot_size_mat[i, j]^(1/3))/2) * min(unit.c(width)),
                                           # height = ((plot_size_mat[i, j]^(1/3))/2) * min(unit.c(height)),
                                           gp = gpar(fill = colors_heatmapbody(plot_data_mat[i, j]), col = "grey70"))
                               }},
                             
                             ## colname parameters
                             show_column_dend = F, cluster_columns = F,
                             top_annotation = colanno_obj, 
                             show_column_names = F,
                             
                             ## row parameters
                             row_gap = unit(x = 0.3, units = "cm"),
                             row_split = row_split_factor, cluster_row_slices = F, cluster_rows = F,
                             left_annotation = rowanno_obj, show_row_names = F,
                             row_title_rot = 0, row_title_gp = gpar(fontsize = 30),
                             show_row_dend = F,
                             show_heatmap_legend = F)
p

# make legend list --------------------------------------------------------
list_lgd = list(
  Legend(col_fun = colors_heatmapbody_legand, 
         title = "Average interaction score", 
         labels_gp = gpar(fontsize=25),
         title_gp = gpar(fontsize = 30),
         legend_width = unit(5, "cm"),
         legend_height = unit(5, "cm"),
         direction = "horizontal"),
  Legend(labels = names(colors_celltype),
         title = "Cell type",
         legend_gp = gpar(fill = colors_celltype), 
         labels_gp = gpar(fontsize=25), title_gp = gpar(fontsize = 30),
         grid_height = unit(1, "cm"),
         grid_width = unit(1, "cm")))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "interactions", ".png")
png(file2write, width =4000 , height = 2600, res = 150)
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()

file2write <- paste0(dir_out, "interactions", ".pdf")
pdf(file2write, width = 25, height = 22)
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()

