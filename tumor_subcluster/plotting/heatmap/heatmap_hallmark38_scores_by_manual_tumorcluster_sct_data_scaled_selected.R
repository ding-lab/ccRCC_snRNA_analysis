# Yige Wu @WashU Jun 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "data.table",
  "stringr",
  "ComplexHeatmap",
  "circlize",
  "RColorBrewer"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    install.packages(pkg_name_tmp, dependencies = T)
  }
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg_name_tmp)
  }
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input meta data
id_metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210423.v1/meta_data.20210423.v1.tsv")
## input pathway scores
scores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/calculate_msigdb_geneset_scores/20210805.v1/MSigDB.Hallmark.tsv")
## input the EMT group later defined
emt_group_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/assign_epithelial_group_bytumorcluster/20211011.v1/Tumorcluster_EpithelialGroup.20211011.v1.tsv")
## input the pathways to plot
genesets_plot_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/pathway/count_ora_sig_genesets_in_up_degs_across_samples/20220606.v1/Count_gene_set_in_up_tumorcluster_degs.20220606.v1.tsv")
## input the annotation for the hallmark gene sets
hallmark_anno_df <- readxl::read_xlsx(path = "./Resources/Knowledge/Databases/MSigDB/Hallmark_gene_sets_summary.xlsx")
## input CNV data
cnv_df <- fread(data.table = F, input = "./Resources/Analysis_Results/copy_number/summarize_cnv_fraction/cnv_fraction_in_tumorcells_per_manualcluster_rm_doublets/20210806.v1/fraction_of_tumorcells_with_cnv_by_gene_by_3state.per_manualsubcluster.20210806.v1.tsv")
## input mutation mapped
driver_mutation_bytumorcluster_df <- fread(data.table = T, input = "./Resources/Analysis_Results/mutation/summarize/count_driver_mutation_mapped_by_intrapatient_tumor_cluster/20220610.v1/Driver_mutation_mapped_per_intrapatienttumorcluster.20220610.v1.tsv")
## input te bulk genomics/methylation events
bulk_sn_omicsprofile_df <- fread(input = "./Resources/Analysis_Results/data_summary/merge_bulk_sn_profiles/20210504.v1/bulk_sn_omics_profile.20210504.v1.tsv", data.table = F)
## input by cluster enrichment assignment
enrich_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/assign_tumorcluster_by_msigdb_geneset_scores/20210805.v1/MsigDB_Hallmark.Top15GeneSets.4Module.Enrichment.tsv")

# preprocess --------------------------------------------------------------
genesets_plot_df <- genesets_plot_df %>%
  mutate(scoregroup_name = paste0(gsub(x = Description, pattern = "HALLMARK_", replacement = ""), "_Score")) %>%
  filter(Description != "SPERMATOGENESIS")
genesets_plot <- genesets_plot_df$scoregroup_name
hallmark_anno_df$`Hallmark Name`[hallmark_anno_df$`Hallmark Name` == "UV_RESPONSE_DOWN"] <- "UV_RESPONSE_DN"
## make cluster id
clustername_df <- data.frame(cluster_name = scores_df$cluster_name)
clustername_df <- clustername_df %>%
  mutate(sampleid = str_split_fixed(string = cluster_name, pattern = "_", n = 2)[,1]) %>%
  mutate(sampleid = gsub(x = sampleid, pattern = "\\.", replacement = "-"))
clustername_df$case <- mapvalues(x = clustername_df$sampleid, from = id_metadata_df$Aliquot.snRNA.WU, to = as.vector(id_metadata_df$Case))
clustername_df <- clustername_df %>%
  filter(case != "C3L-00359") %>%
  filter(!(cluster_name %in% c("C3N.00733.T2_C5", "C3L.01313.T1_C7" , "C3L.01287.T1_C2")))
##
rownames(scores_df) <- scores_df$cluster_name
## prepare CNV data
cnv_df <- cnv_df %>%
  mutate(tumor_subcluster.dataname = gsub(x = tumor_subcluster, pattern = "\\-", replacement = "."))
# format expression data --------------------------------------------------
## get dim names
plot_data_t_mat <- as.matrix(scores_df[,genesets_plot])
plot_data_mat <- t(plot_data_t_mat)
colnames(plot_data_mat) <- scores_df$cluster_name
## make row label
rownames_plot <- rownames(plot_data_mat)
rowlabels_plot <- gsub(x = rownames_plot, pattern = "_Score", replacement = "")
# rowlabels_plot[rowlabels_plot == "EPITHELIAL_MESENCHYMAL_TRANSITION"] <- "EMT"

# make column order -------------------------------------------------------
# col_order_df <- enrich_df %>%
#   arrange(desc(Cell_cycle), desc(Immune), desc(EMT), desc(mTOR))
# plot_data_mat <- plot_data_mat[, col_order_df$cluster_name]
plot_data_mat <- plot_data_mat[, clustername_df$cluster_name]
colnames_plot <- colnames(plot_data_mat)
clusternames_column <- gsub(x = colnames_plot, pattern = "\\.", replacement = "-")
sampleids_column <- str_split_fixed(string = clusternames_column, pattern = "_", n = 2)[,1]

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## colors for gene set scores
summary(as.vector(unlist(plot_data_mat)))
quantile(as.vector(unlist(plot_data_mat)), 0.95)
quantile(as.vector(unlist(plot_data_mat)), 0.05)

colors_heatmapbody = colorRamp2(c(-100, 0, 100), 
                                c("purple", "black", "yellow"))
colors_heatmapbody = colorRamp2(c(quantile(as.vector(unlist(plot_data_mat)), 0.05), 
                                  0, 
                                  quantile(as.vector(unlist(plot_data_mat)), 0.95)), 
                                c("purple", "black", "yellow"))
# colors_heatmapbody = colorRamp2(seq(-100, 100, 20), 
#                                 rev(brewer.pal(n = 11, name = "BrBG")))
# colors for group
# colors_truefalse <- c("black", "grey90")
colors_truefalse <- c("black", "white")
names(colors_truefalse) <- c("TRUE", "FALSE")
color_gridline = "grey90"
colors_enrich_type <- RColorBrewer::brewer.pal(n = 9, name = "Set1")[c(1, 2, 3, 5, 9)]
names(colors_enrich_type) <- paste0(c("EMT", "Cell_cycle", "Immune_signaling", "mTOR", "Other"), "_Module_Enriched")
colors_topbottom <- c("red", "blue", "grey90")
names(colors_topbottom) <- c("top", "bottom", "middle")

# make column annotation --------------------------------------------------
## prepare data
enrich_plot_df <- enrich_df[, c("cluster_name", "Cell_cycle", "Immune", "EMT", "mTOR")]
colnames(enrich_plot_df) <- c("cluster_name", paste0(c("Cell_cycle", "Immune_signaling", "EMT", "mTOR"), "_Module_Enriched"))
tumor_piece_vec <- paste0("T", str_split_fixed(string = colnames_plot, pattern = "\\.T|_", n = 3)[,2])
inflam_score_vec <- scores_df[colnames_plot, "INFLAMMATORY_RESPONSE_Score"]
inflam_assign_vec <- ifelse(inflam_score_vec >= quantile(inflam_score_vec, probs = 0.9), "top",
                            ifelse(inflam_score_vec <= quantile(inflam_score_vec, probs = 0.1), "bottom", "middle"))
# vhl_loss_frac_vec <- mapvalues(x = colnames_plot, from = cnv_df$tumor_subcluster.dataname[cnv_df$gene_symbol == "VHL" & cnv_df$cna_3state == "Loss"], to = as.vector(cnv_df$Fraction[cnv_df$gene_symbol == "VHL" & cnv_df$cna_3state == "Loss"]))
setd2_loss_frac_vec <- mapvalues(x = colnames_plot, from = cnv_df$tumor_subcluster.dataname[cnv_df$gene_symbol == "SETD2" & cnv_df$cna_3state == "Loss"], to = as.vector(cnv_df$Fraction[cnv_df$gene_symbol == "SETD2" & cnv_df$cna_3state == "Loss"]))
setd2_loss_frac_vec[setd2_loss_frac_vec == colnames_plot] <- "0"; setd2_loss_frac_vec <- as.numeric(setd2_loss_frac_vec)
SQSTM1_Gain_frac_vec <- mapvalues(x = colnames_plot, from = cnv_df$tumor_subcluster.dataname[cnv_df$gene_symbol == "SQSTM1" & cnv_df$cna_3state == "Gain"], to = as.vector(cnv_df$Fraction[cnv_df$gene_symbol == "SQSTM1" & cnv_df$cna_3state == "Gain"]))
SQSTM1_Gain_frac_vec[SQSTM1_Gain_frac_vec == colnames_plot] <- "0"; SQSTM1_Gain_frac_vec <- as.numeric(SQSTM1_Gain_frac_vec)
mut_map_vec <- mapvalues(x = clusternames_column, from = driver_mutation_bytumorcluster_df$Cluster_Name, to = as.vector(driver_mutation_bytumorcluster_df$number_cells_w_driver_mutation))
mut_map_vec <- as.character(mut_map_vec != "0")
VHL_bysample_vec <- mapvalues(x = sampleids_column, from = bulk_sn_omicsprofile_df$Aliquot_snRNA_WU, to = as.vector(bulk_sn_omicsprofile_df$Mut.VHL)); VHL_bysample_vec <- as.character(VHL_bysample_vec != "None")
PBRM1_bysample_vec <- mapvalues(x = sampleids_column, from = bulk_sn_omicsprofile_df$Aliquot_snRNA_WU, to = as.vector(bulk_sn_omicsprofile_df$Mut.PBRM1)); PBRM1_bysample_vec <- as.character(PBRM1_bysample_vec != "None")
BAP1_bysample_vec <- mapvalues(x = sampleids_column, from = bulk_sn_omicsprofile_df$Aliquot_snRNA_WU, to = as.vector(bulk_sn_omicsprofile_df$Mut.BAP1)); BAP1_bysample_vec <- as.character(BAP1_bysample_vec != "None")
SETD2_bysample_vec <- mapvalues(x = sampleids_column, from = bulk_sn_omicsprofile_df$Aliquot_snRNA_WU, to = as.vector(bulk_sn_omicsprofile_df$Mut.SETD2)); SETD2_bysample_vec <- as.character(SETD2_bysample_vec != "None")
colanno_obj <- HeatmapAnnotation(chr3p_SETD2_loss_bycluster = anno_simple(x = setd2_loss_frac_vec, col = colorRamp2(seq(0, 1, 0.2), 
                                                                                                          c("white", brewer.pal(n = 6, name = "Blues")[-1]))),
                                 chr5q_SQSTM1_gain_bycluster = anno_simple(x = SQSTM1_Gain_frac_vec, col = colorRamp2(seq(0, 1, 0.2), 
                                                                                                          c("white", brewer.pal(n = 6, name = "Reds")[-1]))),
                                 driver_mutation_bycluster = anno_simple(x = mut_map_vec, col = colors_truefalse[mut_map_vec]),
                                 VHL_mutated_bysample = anno_simple(x = VHL_bysample_vec, col = colors_truefalse[VHL_bysample_vec]),
                                 PBRM1_mutated_bysample = anno_simple(x = PBRM1_bysample_vec, col = colors_truefalse[PBRM1_bysample_vec]),
                                 BAP1_mutated_bysample = anno_simple(x = BAP1_bysample_vec, col = colors_truefalse[BAP1_bysample_vec]),
                                 SETD2_mutated_bysample = anno_simple(x = SETD2_bysample_vec, col = colors_truefalse[SETD2_bysample_vec]),
                                 Inflammatory = anno_simple(x = inflam_assign_vec, col = colors_topbottom[inflam_assign_vec]), annotation_name_side = "left")
# ## merge data
# colanno_df <- enrich_plot_df %>%
#   mutate(EMT_Module_Enriched = (cluster_name %in% emt_group_df$cluster_name[emt_group_df$epithelial_group == "EMT"]))
# rownames(colanno_df) <- colanno_df$cluster_name
# colanno_df <- colanno_df[colnames_plot,-1]
# ## make colors
# colors_scores_list <- list()
# for (colname_tmp in colnames(enrich_plot_df)[-1]) {
#   colanno_df[, colname_tmp] <- as.character(colanno_df[, colname_tmp])
#   colors_tmp <- c(colors_isenriched["FALSE"], colors_enrich_type[colname_tmp])
#   names(colors_tmp) <- c("FALSE", "TRUE")
#   colors_scores_list[[colname_tmp]] <- colors_tmp
# }
# colanno_obj = HeatmapAnnotation(df = colanno_df, col = colors_scores_list,
#                                 annotation_name_gp = gpar(fontsize = 14), annotation_name_side = "left", show_legend = F)


# make row annotation -----------------------------------------------------
freq_de_vec <- mapvalues(x = rownames_plot, from = genesets_plot_df$scoregroup_name, to = as.vector(genesets_plot_df$Freq)); freq_de_vec <- as.numeric(freq_de_vec)
rowanno_obj <- rowAnnotation(Freq_of_DE = anno_barplot(freq_de_vec), annotation_name_side = "top")


# make gene set split -----------------------------------------------------
row_split_vec <- mapvalues(x = rowlabels_plot, from = hallmark_anno_df$`Hallmark Name`, to = as.vector(hallmark_anno_df$`Process Category`))


# make column split -------------------------------------------------------
column_split_vec <- mapvalues(x = colnames_plot, from = clustername_df$cluster_name, to = as.vector(clustername_df$case))
clustercount_df <- clustername_df %>%
  group_by(case) %>%
  summarise(number_clusters = n()) %>%
  arrange(desc(number_clusters))
column_split_factor <- factor(x = column_split_vec, levels = clustercount_df$case)

# plot  ------------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat, 
                             col = colors_heatmapbody,
                             na_col = color_na, #border = "black",
                             ## row
                             show_row_names = T, row_names_gp = gpar(fontsize = 16), row_names_side = "right",
                             show_row_dend = T, row_dend_side = "left", cluster_row_slices = T, 
                             row_split = row_split_vec, row_title_side = "left", row_title_rot = 0, row_title_gp = gpar(fontsize = 16),
                             row_labels = rowlabels_plot, 
                             right_annotation = rowanno_obj,
                             ## column
                             show_column_dend = F, cluster_columns = T, 
                             column_split = column_split_factor, cluster_column_slices = F, column_title_rot = 90,
                             top_annotation = colanno_obj, 
                             show_column_names = F, column_names_side = "top", column_names_gp = gpar(fontsize = 5),
                             show_heatmap_legend = F)
list_lgd = list(
  # Legend(labels = names(colors_isenriched), labels_gp = gpar(fontsize = 14),
  #        title = "Expression\nenriched", title_gp = gpar(fontsize = 14),
  #        legend_gp = gpar(fill = colors_isenriched), border = NA),
  Legend(col_fun = colors_heatmapbody, 
         title = "Gene set score", 
         title_gp = gpar(fontsize = 14),
         labels_gp = gpar(fontsize = 14),
         legend_width = unit(4, "cm"),
         legend_height = unit(4, "cm"),
         direction = "horizontal"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "GeneSetScores", ".pdf")
pdf(file2write, width = 20, height = 13, useDingbats = F)
draw(object = p,
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()


# plot with column name ---------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat, 
                             col = colors_heatmapbody,
                             na_col = color_na, #border = "black",
                             ## row
                             show_row_names = T, row_names_gp = gpar(fontsize = 16), row_names_side = "right",
                             show_row_dend = T, row_dend_side = "left", cluster_row_slices = T, 
                             row_split = row_split_vec, row_title_side = "left", row_title_rot = 0, row_title_gp = gpar(fontsize = 16),
                             row_labels = rowlabels_plot, 
                             right_annotation = rowanno_obj,
                             ## column
                             show_column_dend = F, cluster_columns = T, 
                             column_split = column_split_factor, cluster_column_slices = F, column_title_rot = 90,
                             top_annotation = colanno_obj, 
                             show_column_names = T, column_names_side = "top", column_names_gp = gpar(fontsize = 5),
                             show_heatmap_legend = F)
file2write <- paste0(dir_out, "GeneSetScores", ".png")
png(file2write, width = 2200, height = 1800, res = 150)
draw(object = p,
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()


