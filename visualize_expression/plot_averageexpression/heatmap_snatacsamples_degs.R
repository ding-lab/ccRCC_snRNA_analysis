# Yige Wu @WashU Nov 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
source("https://gist.githubusercontent.com/jokergoo/bfb115200df256eeacb7af302d4e508e/raw/14f315c7418f3458d932ad749850fd515dec413b/word_cloud_grob.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the average expression calculated (SCT)
avgexp_df <- fread(input = "./Resources/Analysis_Results/average_expression/averageexpression_sct_usescale_bycellgroup_w_epithelialcelltypes_byaliquot_on_katmai/20201130.v1/avgexp.SCT.bycellgroup_w_epithelialcelltypes.byaliquot.20201130.v1.tsv", data.table = F)
## input the deg annotation
deg_anno_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/annotate_deg/annotate_deg_by_snatactumorgroup_shared/20201130.v1/DEG_each_tumor_vs_pt.annotated.tsv")
## input idemta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")
## input pathway enrichment
ora_result_df <- fread(data.table = F, input = "./Resources/Analysis_Results/pathway/clustreprofiler/clusterprofiler_wikipathway_deg_snatactumor_vs_normal_alltumor_shared/20201201.v1/ORA_Results.tsv")

# specify parameters ---------------------------------------------------
## specify the samples to show
easyids_snatac <- c("C3N-00733-T1", "C3L-00610-T1", "C3L-01313-T1", "C3L-00416-T2","C3L-00917-T1", "C3L-00088-T1", "C3N-01200-T1", "C3L-00088-T2", "C3L-00448-T1", 
                    # "C3L-01287-T1", "C3L-00079-T1", 
                    "C3L-00088-N", "C3N-01200-N")
aliquots_snatac <- idmetadata_df$Aliquot.snRNA[idmetadata_df$Aliquot.snRNA.WU %in% easyids_snatac]
aliquots_snatac
aliquots_snatac_nat <- idmetadata_df$Aliquot.snRNA[idmetadata_df$Aliquot.snRNA.WU %in% easyids_snatac & idmetadata_df$Sample_Type == "Normal"]
aliquots_snatac_nat
aliquots_snatac_tumor <- idmetadata_df$Aliquot.snRNA[idmetadata_df$Aliquot.snRNA.WU %in% easyids_snatac & idmetadata_df$Sample_Type == "Tumor"]
aliquots_snatac_tumor
## aliquot annotation
sample_anno_df <- data.frame(easyid = c("C3L-00088-N", "C3N-01200-N",
                                        "C3L-00416-T2", "C3L-01313-T1", "C3N-01200-T1",
                                        "C3L-00610-T1", "C3N-00733-T1",
                                        "C3L-00088-T1", "C3L-00088-T2", "C3L-00448-T1", "C3L-00917-T1"),
                             sample_group = c(rep("NAT", 2),
                                              rep("BAP1- tumor", 3),
                                              rep("PBRM1- tumor", 2),
                                              rep("non-mutant tumor", 4)))
## specify genes
deg_filtered_df <- deg_anno_df %>%
  filter(category_byshared %in% c("FALSE_TRUE_FALSE", "FALSE_FALSE_TRUE", "FALSE_TRUE_TRUE", "TRUE_TRUE_TRUE")) %>%
  mutate(category_byshared_label = ifelse(category_byshared == "TRUE_TRUE_TRUE", "all-tumor-shared",
                                          ifelse(category_byshared == "FALSE_TRUE_TRUE", "PBRM-BAP1-shared",
                                                 ifelse(category_byshared == "FALSE_FALSE_TRUE", "PBRM1-shared",
                                                        ifelse(category_byshared == "FALSE_TRUE_FALSE", "BAP1-shared", "Others"))))) %>%
  filter(mean_avg_logFC.bap1mutant > 0.2 | mean_avg_logFC.pbrm1mutant > 0.2) %>%
  arrange(desc(mean_avg_logFC.bap1mutant))
genes2filter <- deg_filtered_df$genesymbol_deg

# format expression data --------------------------------------------------
plot_data_long_df <- avgexp_df %>%
  filter(V1 %in% genes2filter) %>%
  melt() %>%
  mutate(id_bycelltype_byaliquot = gsub(x = variable, pattern = "SCT.", replacement = "")) %>%
  mutate(aliquot = str_split_fixed(string = id_bycelltype_byaliquot, pattern = "_", n = 2)[,1]) %>%
  mutate(celltype_column = str_split_fixed(string = id_bycelltype_byaliquot, pattern = "_", n = 2)[,2])
plot_data_long_df$easyid <- mapvalues(x = plot_data_long_df$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
## filter out non-tumor and NA tumor cluster
plot_data_long_df <- plot_data_long_df %>%
  filter((aliquot %in% aliquots_snatac_nat & celltype_column %in% "Proximal.tubule") | (aliquot %in% aliquots_snatac_tumor & celltype_column == "Tumor.cells")) %>%
  mutate(id_bycluster_byeasyid = paste0(easyid, "_", celltype_column))
## make matrix
plot_data_wide_df <- dcast(data = plot_data_long_df, formula = id_bycluster_byeasyid ~ V1, value.var = "value")
plot_data_mat <- as.matrix(plot_data_wide_df[,-1])
rownames(plot_data_mat) <- plot_data_wide_df$id_bycluster_byeasyid

# get ids -----------------------------------------------------------------
rownames_plot <- rownames(plot_data_mat)
easyids_plot <- str_split_fixed(string = rownames_plot, pattern = "_", n = 2)[,1]
colnames_plot <- colnames(plot_data_mat)

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## make color function for heatmap body colors
color_blue <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[2]
color_red <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
summary(as.vector(plot_data_mat))
colors_heatmapbody = colorRamp2(c(seq(-0.1, -0.5, -0.05),
                                  0,
                                  seq(0.1, 0.5, 0.05)),
                                c(brewer.pal(n = 9, name = "Blues"), "white", brewer.pal(n = 9, name = "YlOrRd")))

# make row split ----------------------------------------------------------
row_split_vec <- mapvalues(x = easyids_plot, from = sample_anno_df$easyid, to = as.vector(sample_anno_df$sample_group))
row_split_factor <- factor(x = row_split_vec, levels = c("BAP1- tumor", "PBRM1- tumor", "non-mutant tumor", "NAT"))

# make column split -------------------------------------------------------
column_split_vec <- mapvalues(x = colnames_plot, from = deg_filtered_df$genesymbol_deg, to = deg_filtered_df$category_byshared_label)
column_split_factor <- factor(x = column_split_vec, levels = c("BAP1-shared", "PBRM1-shared", "PBRM-BAP1-shared", "all-tumor-shared"))

# make column annotation --------------------------------------------------
words_vec <- ora_result_df$Description[ora_result_df$p.adjust < 0.05]
fontsizes_vec <- rep(x = 5, length(words_vec))
gb = word_cloud_grob(text = words_vec, fontsize = fontsizes_vec, max_width = unit(100, "mm"))
gb_h = grobHeight(gb)
gb_w = grobWidth(gb)
# grid.newpage()
# grid.draw(gb)
# grid.rect(width = grobWidth(gb), height = grobHeight(gb), gp = gpar(fill = NA))
# dev.off()
panel_fun = function(index, nm) {
  grid.rect(gp = gpar(fill = "#EEEEEE", col = NA))
  grid.draw(gb)
}
colanno_obj <- HeatmapAnnotation(word_cloud = anno_zoom(align_to = which(column_split_vec == "all-tumor-shared"),
                                                        which = "column", side = "bottom",
                                                        panel_fun = panel_fun, 
                                                        size = gb_h,
                                                        # link_width = gb_w*0.05 + unit(5, "mm"),
                                                        link_width = gb_w*0.1,
                                                        # link_height = gb_w + unit(5, "mm"),
                                                        # link_width = gb_w + unit(5, "mm"),
                                                        link_gp = gpar(fill = "#EEEEEE", col = NA)))

# plot  ------------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat, 
                             col = colors_heatmapbody,
                             na_col = color_na, border = "black",
                             ## row
                             show_row_names = T, row_names_side = "left", row_names_rot = 0, row_title_rot = 0,
                             row_labels = easyids_plot, show_row_dend = F,
                             row_split = row_split_factor, cluster_row_slices = F,
                             # # row_labels = factor_cellgroup,
                             # cluster_row_slices = F, show_row_dend = F, 
                             # ## column
                             show_column_dend = F, cluster_columns = T, 
                             column_split = column_split_factor, show_column_names = F, 
                             # column_title_rot = 15, 
                             column_title_gp = gpar(fontsize = 15),
                             # column_title = NULL,
                             # column_order = column_order_vec,
                             bottom_annotation = colanno_obj,
                             show_heatmap_legend = F)
p

# make legend list --------------------------------------------------------
list_lgd = list(
  Legend(col_fun = colors_heatmapbody, 
         title = "snRNA expression", 
         title_gp = gpar(fontsize = 10),
         labels_gp = gpar(fontsize = 10),
         legend_width = unit(6, "cm"),
         legend_height = unit(4, "cm"),
         direction = "horizontal"))



# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "DEG_expression", ".png")
png(file2write, width = 2000, height = 600, res = 150)
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()
# file2write <- paste0(dir_out, "DEG_expression", ".pdf")
# pdf(file2write, width = 15, height = 5, useDingbats = F)
# draw(object = p, 
#      annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
# dev.off()



