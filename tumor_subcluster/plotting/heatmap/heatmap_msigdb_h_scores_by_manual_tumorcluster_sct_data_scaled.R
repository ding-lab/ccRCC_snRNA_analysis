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
## input meta data
id_metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210423.v1/meta_data.20210423.v1.tsv")
## input pathway scores
scores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/calculate_msigdb_geneset_scores/20210805.v1/MSigDB.Hallmark.tsv")
## input by cluster enrichment assignment
enrich_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/assign_tumorcluster_by_msigdb_geneset_scores/20210805.v1/MsigDB_Hallmark.Top15GeneSets.4Module.Enrichment.tsv")
## input CNV data
cnv_df <- fread(data.table = F, input = "./Resources/Analysis_Results/copy_number/summarize_cnv_fraction/cnv_fraction_in_tumorcells_per_manualcluster_rm_doublets/20210806.v1/fraction_of_tumorcells_with_cnv_by_gene_by_3state.per_manualsubcluster.20210806.v1.tsv")
## input the EMT group later defined
emt_group_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/assign_epithelial_group_bytumorcluster/20211011.v1/Tumorcluster_EpithelialGroup.20211011.v1.tsv")

# preprocess --------------------------------------------------------------
## group gene sets into modules
module1_df <- data.frame(geneset_name = c("HALLMARK_MITOTIC_SPINDLE", "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_DNA_REPAIR", "HALLMARK_MYC_TARGETS_V1"),
                         module_name = "Cell_cycle")
module2_df <- data.frame(geneset_name = c("HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_COMPLEMENT", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_KRAS_SIGNALING_UP"),
                         module_name = "Immune_signaling")
module3_df <- data.frame(geneset_name = c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),
                         module_name = "EMT")
# module4_df <- data.frame(geneset_name = c("HALLMARK_UV_RESPONSE_DN", "HALLMARK_MTORC1_SIGNALING"),
#                          module_name = "mTOR")
module4_df <- data.frame(geneset_name = c("HALLMARK_MTORC1_SIGNALING"),
                         module_name = "mTOR")
modules_df <- rbind(module1_df, module2_df, module3_df, module4_df)
modules_df <- modules_df %>%
  mutate(scoregroup_name =  paste0(gsub(x = geneset_name, pattern = "HALLMARK_", replacement = ""), "_Score"))

# format expression data --------------------------------------------------
## get dim names
scorenames <- c(modules_df$scoregroup_name, "UV_RESPONSE_DN_Score", "HYPOXIA_Score", "TNFA_SIGNALING_VIA_NFKB_Score"); scorenames <- unique(scorenames)
plot_data_t_mat <- as.matrix(scores_df[,scorenames])
plot_data_mat <- t(plot_data_t_mat)
colnames(plot_data_mat) <- scores_df$cluster_name
## make row label
rownames_plot <- rownames(plot_data_mat)
rowlabels_plot <- gsub(x = rownames_plot, pattern = "_Score", replacement = "")
rowlabels_plot[rowlabels_plot == "EPITHELIAL_MESENCHYMAL_TRANSITION"] <- "EMT"


# make column order -------------------------------------------------------
col_order_df <- enrich_df %>%
  arrange(desc(Cell_cycle), desc(Immune), desc(EMT), desc(mTOR))
plot_data_mat <- plot_data_mat[, col_order_df$cluster_name]
colnames_plot <- colnames(plot_data_mat)

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## colors for gene set scores
color_orange <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[5]
color_purple <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[4]
colors_heatmapbody = colorRamp2(c(-100, 0, 100), 
                                c(color_purple, "white", color_orange))
colors_heatmapbody = colorRamp2(seq(-100, 100, 20), 
                                rev(brewer.pal(n = 11, name = "BrBG")))
# colors for group
colors_isenriched <- c("black", "grey90")
names(colors_isenriched) <- c("TRUE", "FALSE")
color_gridline = "grey90"
colors_enrich_type <- RColorBrewer::brewer.pal(n = 9, name = "Set1")[c(1, 2, 3, 5, 9)]
names(colors_enrich_type) <- paste0(c("EMT", "Cell_cycle", "Immune_signaling", "mTOR", "Other"), "_Module_Enriched")

# make column annotation --------------------------------------------------
## prepare data
enrich_plot_df <- enrich_df[, c("cluster_name", "Cell_cycle", "Immune", "EMT", "mTOR")]
colnames(enrich_plot_df) <- c("cluster_name", paste0(c("Cell_cycle", "Immune_signaling", "EMT", "mTOR"), "_Module_Enriched"))

## merge data
colanno_df <- enrich_plot_df %>%
  mutate(EMT_Module_Enriched = (cluster_name %in% emt_group_df$cluster_name[emt_group_df$epithelial_group == "EMT"]))
rownames(colanno_df) <- colanno_df$cluster_name
colanno_df <- colanno_df[colnames_plot,-1]
## make colors
colors_scores_list <- list()
for (colname_tmp in colnames(enrich_plot_df)[-1]) {
  colanno_df[, colname_tmp] <- as.character(colanno_df[, colname_tmp])
  colors_tmp <- c(colors_isenriched["FALSE"], colors_enrich_type[colname_tmp])
  names(colors_tmp) <- c("FALSE", "TRUE")
  colors_scores_list[[colname_tmp]] <- colors_tmp
}
colanno_obj = HeatmapAnnotation(df = colanno_df, col = colors_scores_list,
                                annotation_name_gp = gpar(fontsize = 14), annotation_name_side = "left", show_legend = F)

# plot  ------------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat, 
                             col = colors_heatmapbody,
                             na_col = color_na, #border = "black",
                             ## row
                             show_row_names = T, row_names_gp = gpar(fontsize = 14), row_names_side = "left",
                             show_row_dend = F, cluster_row_slices = T, row_labels = rowlabels_plot,
                             ## column
                             show_column_dend = F, cluster_columns = F,
                             top_annotation = colanno_obj,
                             show_column_names = F, column_names_side = "top", column_names_gp = gpar(fontsize = 5),
                             show_heatmap_legend = F)
list_lgd = list(
  Legend(labels = names(colors_isenriched), labels_gp = gpar(fontsize = 14),
         title = "Expression\nenriched", title_gp = gpar(fontsize = 14),
         legend_gp = gpar(fill = colors_isenriched), border = NA),
  Legend(col_fun = colors_heatmapbody, 
         title = "Gene set score", 
         title_gp = gpar(fontsize = 14),
         labels_gp = gpar(fontsize = 14),
         legend_width = unit(4, "cm"),
         legend_height = unit(4, "cm"),
         direction = "vertical"))

file2write <- paste0(dir_out, "GeneSetScores", ".pdf")
pdf(file2write, width = 10, height = 4, useDingbats = F)
draw(object = p,
     annotation_legend_side = "right", annotation_legend_list = list_lgd)
dev.off()


# plot with column name ---------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat, 
                             col = colors_heatmapbody,
                             na_col = color_na, #border = "black",
                             ## row
                             show_row_names = T, row_names_gp = gpar(fontsize = 14), row_names_side = "left",
                             show_row_dend = F, cluster_row_slices = T, row_labels = rowlabels_plot,
                             ## column
                             show_column_dend = F, cluster_columns = F,
                             top_annotation = colanno_obj,
                             show_column_names = T, column_names_side = "top", column_names_gp = gpar(fontsize = 5),
                             show_heatmap_legend = F)
file2write <- paste0(dir_out, "GeneSetScores", ".png")
png(file2write, width = 1300, height = 900, res = 150)
draw(object = p, 
     annotation_legend_side = "top", annotation_legend_list = list_lgd)
dev.off()


