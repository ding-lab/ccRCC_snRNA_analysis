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
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the average expression calculated (SCT)
avgexp_df <- fread(input = "./Resources/Analysis_Results/average_expression/averageexpression_sct_usescale_bycellgroup_w_epithelialcelltypes_byaliquot_on_katmai/20201130.v1/avgexp.SCT.bycellgroup_w_epithelialcelltypes.byaliquot.20201130.v1.tsv", data.table = F)
## input the deg annotation
deg_anno_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/pbrm1_bap1_vs_non_mutants/summarize_degs/summarize_PBRM1_BAP1_DEGs/20210405.v1/BAP1_PBRM1_DEGs.Num_samples.20210405.v1.tsv")
## input idemta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")
## input sample annotaiton
sample_anno_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/mutation/annotate_cptac_sample_by_pbrm1_bap1_mutation/20210310.v1/PBRM1_BAP1_Mutation_Status_By_Case.20210310.v1.tsv")

# specify parameters ---------------------------------------------------
aliquots_nat <- idmetadata_df$Aliquot.snRNA[idmetadata_df$Sample_Type == "Normal"]
aliquots_nat
aliquots_tumor <- idmetadata_df$Aliquot.snRNA[idmetadata_df$Sample_Type == "Tumor"]
aliquots_tumor

## specify genes
deg_filtered_df <- deg_anno_df %>%
  filter(PBRM1_vs_NonMutants_snRNA != "Inconsistent")
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
  filter((aliquot %in% aliquots_nat & celltype_column %in% "Proximal.tubule") | (aliquot %in% aliquots_tumor & celltype_column == "Tumor.cells")) %>%
  mutate(id_bycluster_byeasyid = paste0(easyid, "_", celltype_column))
## make matrix
plot_data_wide_df <- dcast(data = plot_data_long_df, formula = id_bycluster_byeasyid ~ V1, value.var = "value")
plot_data_mat <- as.matrix(plot_data_wide_df[,-1])
rownames(plot_data_mat) <- plot_data_wide_df$id_bycluster_byeasyid

# get ids -----------------------------------------------------------------
rownames_plot <- rownames(plot_data_mat)
easyids_plot <- str_split_fixed(string = rownames_plot, pattern = "_", n = 2)[,1]
colnames_plot <- colnames(plot_data_mat)
caseids_plot <- mapvalues(x = easyids_plot, from = idmetadata_df$Aliquot.snRNA.WU, to = as.vector(idmetadata_df$Case))

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
row_split_vec <- mapvalues(x = caseids_plot, from = sample_anno_df$Case, to = as.vector(sample_anno_df$mutation_category_sim))
row_split_vec[grepl(x = easyids_plot, pattern = "\\-N")] <- "NAT"
row_split_factor <- factor(x = row_split_vec, levels = c("PBRM1 mutated", "BAP1 mutated", "Non-mutants", "Both mutated", "NAT"))

# make column split -------------------------------------------------------
column_split_vec <- mapvalues(x = colnames_plot, from = deg_filtered_df$genesymbol_deg, to = deg_filtered_df$PBRM1_vs_NonMutants_snRNA)

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
                             column_split = column_split_vec, show_column_names = F, 
                             # column_title_rot = 15, 
                             column_title_gp = gpar(fontsize = 15),
                             # column_title = NULL,
                             # column_order = column_order_vec,
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
file2write <- paste0(dir_out, "PBRM1_vs_NonMutant.DEG_expression", ".png")
png(file2write, width = 2000, height = 800, res = 150)
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()



