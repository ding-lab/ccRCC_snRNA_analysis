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
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input summary for cell-cell interactin
summary_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/summarize_interactions/rank_across_celltypes_by_pair_by_avgsigmean/20200925.v3/cellphonedb.summary_across_celltypes_by_pair.min5.tsv")

# correct -----------------------------------------------------------------
summary_df$celltypes.source2target[summary_df$interacting_pair == "HLA-G_LILRB2"] <- "Tumor cells->Macrophages"
summary_df$gene.source[summary_df$interacting_pair == "HLA-G_LILRB2"] <- "HLA-G"
summary_df$gene.target[summary_df$interacting_pair == "HLA-G_LILRB2"] <- "LILRB2"
summary_df$Cell_group.source[summary_df$interacting_pair == "HLA-G_LILRB2"] <- "Tumor cells"


# specify pairs to filter -------------------------------------------------
paired_cellgroups.general_process <- "Tumor&Immune"
summary_filtered_df <- summary_df %>%
  filter(paired_cellgroups.general == paired_cellgroups.general_process) %>%
  filter(rank_byinteraction_acrosscelltypes == 1) %>%
  filter(!(Cell_type.source %in% c("Mixed myeloid/lymphoid")) & !(Cell_type.target %in% c("Mixed myeloid/lymphoid"))) %>%
  filter(!(interacting_pair %in% c("AXL_IL15RA"))) %>%
  mutate(paired_celltypes_group = paste0("Tumor&",ifelse(Cell_type.source == "Tumor cells", Cell_type.target, Cell_type.source))) %>%
  arrange(desc(avg_sig_mean), number_sig_cases)
interacting_pairs_process <- summary_filtered_df$interacting_pair
interacting_pairs_process
paired_celltypes_process <- summary_filtered_df$paired_celltypes
length(interacting_pairs_process)

# make data matrix for the heatmap body -------------------------------------------
plot_data_df <- summary_df %>%
  filter(interacting_pair %in% interacting_pairs_process) %>%
  filter(paired_celltypes %in% paired_celltypes_process) %>%
  filter(paired_cellgroups.general == paired_cellgroups.general_process)
plot_data_wide_df <- dcast(data = plot_data_df, formula = interacting_pair ~ celltypes.source2target, value.var = "avg_sig_mean")
plot_data_mat <- as.matrix(plot_data_wide_df[, -1])
rownames(plot_data_mat) <- plot_data_wide_df$interacting_pair
plot_data_mat <- plot_data_mat[interacting_pairs_process,]
## make non-AN values that are not 1 as 0 as
# plot_data_mat[!is.na(plot_data_mat) & plot_data_mat != 1] <- 0
plot_data_mat[is.na(plot_data_mat)] <- -1

# get ids -----------------------------------------------------------------
interacting_pair_plot <- rownames(plot_data_mat)
celltypes.source2target_plot <- colnames(plot_data_mat)

# specify colors ----------------------------------------------------------
color_red <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
color_orange <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[5]
color_yellow <- RColorBrewer::brewer.pal(n = 9, name = "YlGnBu")[1]
## specify color for NA values
color_na <- "grey50"
## make color function for heatmap body colors
colors_heatmapbody <- c(color_na, "white", color_red)
names(colors_heatmapbody) <- c(-1, 0, 1)
colors_heatmapbody
# colors_heatmapbody = colorRamp2(c(-1, 
#                                   0, 
#                                   1), 
#                                 c("black", "white", color_red))
## cell type colors
celltypes_uniq <- unique(c(unique(summary_filtered_df$Cell_type.target), unique(summary_filtered_df$Cell_type.source)))
celltypes_uniq
RColorBrewer::display.brewer.all()
colors_celltype <- c(cellgroup_colors["Tumor cells"], RColorBrewer::brewer.pal(n = (length(celltypes_uniq)-1), name = "Set3"))
names(colors_celltype) <- c("Tumor cells", celltypes_uniq[celltypes_uniq != "Tumor cells"])
colors_celltype
## make colors for average sig means
colors_heatmapbody <- colorRamp2(c(-1, 0, 2, 4, 6, 8, 10, 12, 14, 16), 
                                 c(color_na, brewer.pal(n = 9, name = "YlOrRd")))

# make row annotation -----------------------------------------------------
avg_sig_mean_vec <- mapvalues(x = interacting_pair_plot, from = summary_filtered_df$interacting_pair, to = as.vector(summary_filtered_df$avg_sig_mean))
avg_sig_mean_vec <- as.numeric(avg_sig_mean_vec)
avg_sig_mean_vec
summary(avg_sig_mean_vec)
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
## make interaction texts
interaction_direction_vec <- mapvalues(x = interacting_pair_plot, from = summary_filtered_df$interacting_pair, to = as.vector(summary_filtered_df$is_ligand_receptor))
interaction_direction_text <- ifelse(interaction_direction_vec == "TRUE", "->", "&")
## make interaction direciton
# make annotation object
rowanno_obj1 <- rowAnnotation(
  # Celltype_ligand = anno_simple(x = vector(mode = "numeric", length = length(celltype_source_vec)), 
  #                               gp = gpar(fill = "white"),
  #                               pch = 21, pt_gp = gpar(fill = colors_celltype[celltype_source_vec])),
  Gene_source = anno_text(x = gene_source_vec, gp = gpar(color = "black")),
  # Celltype_receptor = anno_simple(x = vector(mode = "numeric", length = length(celltype_target_vec)), 
  #                                 gp = gpar(fill = "white"),
  #                                 pch = 21, pt_gp = gpar(fill = colors_celltype[celltype_target_vec])),
  Gene_target = anno_text(x = gene_target_vec),
  annotation_name_side = "top")
rowanno_obj2 <- rowAnnotation(Number_samples_significant = anno_barplot(x = number_sig_samples_vec, width = unit(x = 3, "cm")),
                              annotation_name_side = "top")

# make column annotation -----------------------------------------------------
col_celltypes.source_vec <- str_split_fixed(string = celltypes.source2target_plot, pattern = '->|&', n = 2)[,1]
col_celltypes.target_vec <- str_split_fixed(string = celltypes.source2target_plot, pattern = '->|&', n = 2)[,2]
## get the non-tumor cell types
nontumor_celltypes_col_vec <- col_celltypes.target_vec
nontumor_celltypes_col_vec[nontumor_celltypes_col_vec == "Tumor cells"] <- col_celltypes.source_vec[nontumor_celltypes_col_vec == "Tumor cells"]
## get direction
direction_pch_vec <- ifelse(col_celltypes.source_vec == "Tumor cells", 25, 24)
## make column annotation object
colanno_obj <- HeatmapAnnotation(Celltype1 = anno_simple(x = vector(mode = "numeric", length = length(col_celltypes.source_vec)), 
                                                         gp = gpar(fill = "white"), height = unit(6, "mm"), 
                                                         pch = 21, pt_size = unit(6, "mm"), pt_gp = gpar(fill = colors_celltype["Tumor cells"])), 
                                 Direction = anno_simple(x = vector(mode = "numeric", length = length(col_celltypes.source_vec)), 
                                                         gp = gpar(fill = "white"), height = unit(5, "mm"), 
                                                         pch = direction_pch_vec, pt_size = unit(5, "mm"), pt_gp = gpar(fill = "black")),
                                 Celltype2 = anno_simple(x = vector(mode = "numeric", length = length(col_celltypes.source_vec)), 
                                                         gp = gpar(fill = "white"), height = unit(6, "mm"), 
                                                         pch = 21, pt_size = unit(6, "mm"), pt_gp = gpar(fill = colors_celltype[nontumor_celltypes_col_vec])),
                                 annotation_name_side = "left")


# make column split -------------------------------------------------------
nontumor_celltypes_col_vec
col_split_vector <- mapvalues(x = nontumor_celltypes_col_vec, 
                              from = c("Macrophages", "Macrophages proliferating", "TRM",
                                       "B-cells", "Plasma",
                                       "CD4 CTL",
                                       "NK cells strong", "NK cells weak",
                                       "cDC",
                                       "CD4/CD8 proliferating"),
                              to = c("Macrophages", "Macrophages", "Macrophages", 
                                     "B-cells", "B-cells",
                                     "CD4 T-cells",
                                     "NK cells", "NK cells",
                                     "DC",
                                     "CD4/CD8 proliferating"))
col_split_vector

# plot hetamp body --------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat,
                             col = colors_heatmapbody,
                             border = "black", 
                             ## row parameters
                             show_row_dend = F, cluster_rows = F,
                             left_annotation = rowanno_obj1, right_annotation = rowanno_obj2,
                             show_row_names = F,
                             ## column parameters
                             column_split = col_split_vector,
                             top_annotation = colanno_obj, show_column_names = F,
                             column_names_side = "top", column_title_rot = 90,
                             # column_title = NULL,
                             show_column_dend = F,
                             show_heatmap_legend = F)
p

# make legend list --------------------------------------------------------
list_lgd = list(
  Legend(col_fun = colors_heatmapbody, 
         title = "Interaction strength score", 
         labels_gp = gpar(fontsize=12), title_gp = gpar(fontsize = 15),
         legend_width = unit(4, "cm"),
         legend_height = unit(3, "cm"),
         direction = "horizontal"),
  Legend(labels = names(colors_celltype),
         title = "Cell_type",
         legend_gp = gpar(fill = colors_celltype), 
         labels_gp = gpar(fontsize=12), title_gp = gpar(fontsize = 15),
         grid_height = unit(1, "cm"),
         grid_width = unit(1, "cm")))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "immune_tumor_interactions.", ".png")
png(file2write, width = 1700, height = 1600, res = 150)
draw(object = p, 
     annotation_legend_side = "left", annotation_legend_list = list_lgd)
dev.off()

file2write <- paste0(dir_out, "immune_tumor_interactions.", ".pdf")
pdf(file2write, width = 15, height = 20)
draw(object = p, 
     annotation_legend_side = "left", annotation_legend_list = list_lgd)
dev.off()

