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

# specify pairs to filter -------------------------------------------------
paired_cellgroups.general_process <- "Tumor&Immune"
summary_filtered_df <- summary_df %>%
  filter(paired_cellgroups.general == paired_cellgroups.general_process) %>%
  filter(rank_byinteraction_acrosscelltypes == 1) %>%
  mutate(paired_celltypes_group = paste0("Tumor&",ifelse(Cell_type.source == "Tumor cells", Cell_type.target, Cell_type.source))) %>%
  arrange(paired_celltypes_group, desc(paired_cellgroups.detailed), desc(avg_sig_mean), number_sig_cases)
interacting_pairs_process <- summary_filtered_df$interacting_pair
interacting_pairs_process
length(interacting_pairs_process)
# make data matrix for the heatmap body -------------------------------------------
plot_data_df <- summary_df %>%
  filter(interacting_pair %in% interacting_pairs_process) %>%
  filter(paired_cellgroups.general == paired_cellgroups.general_process)
plot_data_wide_df <- dcast(data = plot_data_df, formula = interacting_pair ~ celltypes.source2target, value.var = "rank_byinteraction_acrosscelltypes")
plot_data_mat <- as.matrix(plot_data_wide_df[, -1])
rownames(plot_data_mat) <- plot_data_wide_df$interacting_pair
plot_data_mat <- plot_data_mat[interacting_pairs_process,]
## make non-AN values that are not 1 as 0 as
plot_data_mat[!is.na(plot_data_mat) & plot_data_mat != 1] <- 0
plot_data_mat[is.na(plot_data_mat)] <- -1

# get ids -----------------------------------------------------------------
interacting_pair_plot <- rownames(plot_data_mat)
celltypes.source2target_plot <- colnames(plot_data_mat)
paired_cellgroups_detailed <- mapvalues(x = celltypes.source2target_plot, from = summary_filtered_df$celltypes.source2target, to = as.vector(summary_filtered_df$paired_celltypes_group)) 
paired_cellgroups_detailed

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
colors_sig_means <- colorRamp2(c(0, 2, 4, 6, 8, 10), 
                               brewer.pal(n = 6, name = "YlGn"))


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
## map pathways
pathway_vec <- vector(mode = "character", length = length(gene_source_vec))
pathway_vec[gene_source_vec %in% c("VEGFA", "VEGFB", "VEGFC") | gene_target_vec %in% c("FLT1", "KDR", "NRP1", "NRP2")] <- "VEGFR"
pathway_vec[gene_target_vec %in% c("EGFR")] <- "EGFR"
pathway_vec[gene_target_vec %in% c("ERBB3", "ERBB2")] <- "ERBB"
pathway_vec[gene_target_vec %in% c("EPHA2", "EPHA4", "EPHA7")] <- "Ephrin"
pathway_vec[gene_target_vec %in% c("FGFR1", "FGFR2", "FGFR3", "FGFR4")] <- "FGFR"
pathway_vec[gene_target_vec %in% c("LGR4")] <- "WNT"
pathway_vec[gene_source_vec %in% c("LIF", "FAM3C")] <- "EMT"
pathway_vec[gene_source_vec %in% c("BMP6")] <- "TGFb"
pathway_vec[gene_target_vec %in% summary_filtered_df$gene.target[summary_filtered_df$is_integrin]] <- "Integrin"
## make interaction direciton
# make annotation object
rowanno_obj1 <- rowAnnotation(Gene_source = anno_text(x = gene_source_vec, gp = gpar(fill = colors_celltype[celltype_source_vec])),
                              Gene_target = anno_text(x = paste0(interaction_direction_text,gene_target_vec), gp = gpar(fill = colors_celltype[celltype_target_vec])),
                              Pathway = anno_text(x = pathway_vec),
                              annotation_name_side = "top")
rowanno_obj2 <- rowAnnotation(Interaction_strength_top = anno_simple(x = avg_sig_mean_vec, col = colors_sig_means),
                              Number_samples_significant = anno_barplot(x = number_sig_samples_vec, width = unit(x = 3, "cm")),
                              annotation_name_side = "top")
# make column annotation -----------------------------------------------------
paired_cellgroups_detailed <- mapvalues(x = celltypes.source2target_plot, from = summary_filtered_df$celltypes.source2target, to = as.vector(summary_filtered_df$paired_celltypes_group)) 
col_celltypes.source_vec <- mapvalues(x = celltypes.source2target_plot, from = summary_filtered_df$celltypes.source2target, to = as.vector(summary_filtered_df$Cell_type.source))
col_celltypes.target_vec <- mapvalues(x = celltypes.source2target_plot, from = summary_filtered_df$celltypes.source2target, to = as.vector(summary_filtered_df$Cell_type.target))
colanno_obj <- HeatmapAnnotation(Celltype_source = anno_simple(x = paste0("<-", col_celltypes.source_vec), gp = gpar(fill = colors_celltype[col_celltypes.source_vec])),
                                 Celltype_target = anno_simple(x = col_celltypes.target_vec, gp = gpar(fill = colors_celltype[col_celltypes.target_vec])))

# plot hetamp body --------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat,
                             col = colors_heatmapbody,
                             border = "black", 
                             ## row parameters
                             # row_names_side = "left",
                             # row_names_max_width = unit(10, "cm"),
                             show_row_dend = F, cluster_rows = F,
                             left_annotation = rowanno_obj1, right_annotation = rowanno_obj2,
                             show_row_names = F,
                             ## column parameters
                             column_split = paired_cellgroups_detailed,
                             top_annotation = colanno_obj, show_column_names = F,
                             column_names_side = "top", column_title = NULL,
                             show_column_dend = F,
                             show_heatmap_legend = F)
p

# make legend list --------------------------------------------------------
list_lgd = list(
  Legend(labels = c("significant over 5 samples, rank #1",
                    "significant over 5 samples, not rank #1",
                    "Not significant over 5 samples"), 
         legend_gp = gpar(fill = colors_heatmapbody[c("1", "0", "-1")]), 
         labels_gp = gpar(fontsize=12), title_gp = gpar(fontsize = 15),
         title = "Rank of interaction strength\nacross samples across cell types", 
         grid_height = unit(1, "cm"),
         grid_width = unit(1, "cm"),
         direction = "horizontal"),
  Legend(col_fun = colors_sig_means, 
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
file2write <- paste0(dir_out, "stroma_tumor_interactions.", ".png")
png(file2write, width = 1500, height = 2000, res = 150)
draw(object = p, 
     annotation_legend_side = "left", annotation_legend_list = list_lgd)
dev.off()

file2write <- paste0(dir_out, "stroma_tumor_interactions.", ".pdf")
pdf(file2write, width = 15, height = 20)
draw(object = p, 
     annotation_legend_side = "left", annotation_legend_list = list_lgd)
dev.off()

