# Yige Wu @WashU Nov 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the average expression calculated (SCT)
avgexp_long_df <- fread(input = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/Enriched_Motifs/Tumor_vs_NormalPT/Motif_score_perCell_group.tsv", data.table = F)
## input the dam annotation
dam_selected_df1 <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/enriched_motifs/summarize_tf_motifs/summarize_sig_up_tf_diff_across_snatac_tumors/20201223.v2/DAMs_top.20201223.v2.tsv")
dam_selected_df2 <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/enriched_motifs/summarize_tf_motifs/summarize_sig_down_tf_diff_across_snatac_tumors/20201227.v1/DAMs_top_down.20201227.v1.tsv")
## input idemta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")

# specify parameters ---------------------------------------------------
## specify the samples to show
easyids_snatac <- c("C3N-00733-T1", "C3L-00610-T1", "C3L-01313-T1", "C3L-00416-T2","C3L-00917-T1", "C3L-00088-T1", "C3N-01200-T1", "C3L-00088-T2", 
                    "C3L-00448-T1",
                    "C3L-01287-T1", "C3L-00079-T1",
                    "C3L-00088-N", "C3N-01200-N")
aliquots_snatac <- idmetadata_df$Aliquot.snRNA[idmetadata_df$Aliquot.snRNA.WU %in% easyids_snatac]
aliquots_snatac
aliquots_snatac_nat <- idmetadata_df$Aliquot.snRNA[idmetadata_df$Aliquot.snRNA.WU %in% easyids_snatac & idmetadata_df$Sample_Type == "Normal"]
aliquots_snatac_nat
aliquots_snatac_tumor <- idmetadata_df$Aliquot.snRNA[idmetadata_df$Aliquot.snRNA.WU %in% easyids_snatac & idmetadata_df$Sample_Type == "Tumor"]
aliquots_snatac_tumor
## aliquot annotation
sample_anno_df <- data.frame(easyid = c("C3L-00088-N", "C3N-01200-N",
                                        "C3L-01313-T1", "C3N-01200-T1", "C3L-01287-T1", "C3L-00416-T2", 
                                        "C3L-00610-T1", "C3N-00733-T1", "C3L-00079-T1", "C3L-00416-T2", 
                                        "C3L-00088-T1", "C3L-00088-T2", "C3L-00448-T1", "C3L-00917-T1"),
                             sample_group = c(rep("NAT", 2),
                                              rep("BAP1-mutant tumor", 4),
                                              rep("PBRM1-mutant tumor", 4),
                                              rep("non-mutant tumor", 4)))
## specify genes
dam_anno_df <- rbind(dam_selected_df1 %>%
                           filter(TF_category == "tumor-common-up") %>%
                           select(TF_Name, TF_category),
                         dam_selected_df2 %>%
                           select(TF_Name, TF_category))
genes2filter <- dam_anno_df$TF_Name

# format expression data --------------------------------------------------
plot_data_long_df <- avgexp_long_df %>%
  filter(TF_Name %in% genes2filter) %>%
  filter(cell_type %in% easyids_snatac)
## filter out non-tumor and NA tumor cluster
## make matrix
plot_data_wide_df <- dcast(data = plot_data_long_df, formula = cell_type ~ TF_Name, value.var = "mean_score")
plot_data_mat <- as.matrix(plot_data_wide_df[,-1])
rownames(plot_data_mat) <- plot_data_wide_df$cell_type

# get ids -----------------------------------------------------------------
rownames_plot <- rownames(plot_data_mat)
easyids_plot <- rownames_plot
colnames_plot <- colnames(plot_data_mat)

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## make color function for heatmap body colors
color_blue <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[2]
color_red <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
summary(as.vector(plot_data_mat))
colors_heatmapbody = colorRamp2(c(rev(seq(-0.1, -0.5, -0.05)),
                                  0,
                                  c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1, 2, 3)),
                                c(rev(brewer.pal(n = 9, name = "Blues")), "white", brewer.pal(n = 9, name = "YlOrRd")))
## make colors for direction
colors_direction <- c("up" = color_red, "down" = color_blue)

# make row split ----------------------------------------------------------
row_split_vec <- mapvalues(x = easyids_plot, from = sample_anno_df$easyid, to = as.vector(sample_anno_df$sample_group))
row_split_factor <- factor(x = row_split_vec, levels = c("non-mutant tumor","BAP1-mutant tumor", "PBRM1-mutant tumor",  "NAT"))

# make column split -------------------------------------------------------
column_split_vec <- mapvalues(x = colnames_plot, from = dam_anno_df$TF_Name, to = as.vector(dam_anno_df$TF_category))
column_split_factor <- factor(x = column_split_vec, levels = c("tumor-common-up", "tumor-common-down"))

# make column annotation --------------------------------------------------
## direction
dam_direction_vec <- mapvalues(x = colnames_plot, from = dam_anno_df$TF_Name, to = dam_anno_df$direction_shared)
## filter motifs to show
dam_anno_df$mean_diff_acrosstumorgroups <- rowMeans(x = dam_anno_df[, c("mean_diff.nonmutant","mean_diff.bap1mutant","mean_diff.pbrm1mutant")], na.rm = T)
dam_top_df <- dam_anno_df %>%
  mutate(abs_mean_diff = abs(mean_diff_acrosstumorgroups)) %>%
  group_by(category_byshared_label, direction_shared) %>%
  top_n(n = 10, wt = abs_mean_diff)
index_highlight <- which(colnames_plot %in% dam_top_df$TF_Name)
texts_highlight <- colnames_plot[index_highlight]
index_upmotifs <- which(colnames_plot %in% dam_top_df$TF_Name[dam_top_df$direction_shared == "up"])
texts_upmotifs <- colnames_plot[index_upmotifs]
index_downmotifs <- which(colnames_plot %in% dam_top_df$TF_Name[dam_top_df$direction_shared == "down"])
texts_downmotifs <- colnames_plot[index_downmotifs]
colanno_obj1 <- HeatmapAnnotation(Motifs_up = anno_mark(at = index_upmotifs, labels = texts_upmotifs, labels_gp = gpar(fontsize = 8, color = "red"), side = "top"),
                                  Tumor_vs_PT = anno_simple(x = dam_direction_vec, col = colors_direction, height = unit(0.2, "cm")),
                                  annotation_name_side = "left")
colanno_obj2 <- HeatmapAnnotation(Motifs_down = anno_mark(at = index_downmotifs, labels = texts_downmotifs, labels_gp = gpar(fontsize = 8), side = "bottom"),
                                  annotation_name_side = "left")

# plot  ------------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat, 
                             col = colors_heatmapbody,
                             na_col = color_na,
                             ## row
                             show_row_names = T, row_names_side = "left", row_names_rot = 180, row_title_rot = 0,
                             row_labels = easyids_plot, show_row_dend = F, row_names_gp = gpar(fontsize = 15),
                             row_split = row_split_factor, cluster_row_slices = F, 
                             # # row_labels = factor_cellgroup,
                             # cluster_row_slices = F, show_row_dend = F, 
                             # ## column
                             show_column_dend = F, cluster_columns = T, 
                             column_split = column_split_factor,
                             show_column_names = T,  column_names_side = "top", column_names_gp = gpar(fontsize = 17),
                             # column_title_rot = 15, 
                             column_title_gp = gpar(fontsize = 15), column_title_rot = 0, column_title_side = "bottom",
                             # column_title = NULL,
                             # column_order = column_order_vec,
                             show_heatmap_legend = F)
p

# make legend list --------------------------------------------------------
list_lgd = list(
  Legend(col_fun = colors_heatmapbody, 
         title = "Motif score", 
         title_gp = gpar(fontsize = 10),
         labels_gp = gpar(fontsize = 10),
         legend_width = unit(3, "cm"),
         legend_height = unit(4, "cm"),
         direction = "horizontal"))



# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "dam_expression", ".png")
png(file2write, width = 1500, height = 1000, res = 150)
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()
file2write <- paste0(dir_out, "dam_expression", ".pdf")
pdf(file2write, width = 9, height = 5.5, useDingbats = F)
draw(object = p,
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()



