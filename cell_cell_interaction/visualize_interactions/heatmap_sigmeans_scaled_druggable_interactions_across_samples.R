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
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input summary for cell-cell interactin
summary_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/filter_interactions/filter_druggable_related_pair_celltypes/20200924.v1/filtered_druggale_related_pair_celltypes.tsv")
## input cell-cell interaction values by case
cellphone_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/other/scale_pair_cell_types_vs_sample_sig_mean/20200924.v1/cellphonedb.pair_cell.types_vs_sample.sig_mean.scaled_by_sample_bypair..tsv")
## input sample meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")

# specify pairs to filter -------------------------------------------------
nrow(summary_df)
summary_filtered_df <- summary_df %>%
  filter(therapy_category != "Other") %>%
  filter(!(interacting_pair %in% c("FLT1 complex_VEGFA", "FLT1 complex_VEGFB")))
pair_cell.types_process <- summary_filtered_df$interaction_celltypes
pair_cell.types_process

# make data matrix for the heatmap body -------------------------------------------
plot_data_df <- cellphone_df %>%
  filter(V1 %in% pair_cell.types_process)
plot_data_mat <- as.matrix(plot_data_df[, -1])
rownames(plot_data_mat) <- plot_data_df$V1
plot_data_mat[1:5, 1:10]
plot_data_mat[is.na(plot_data_mat)] <- -2
ids_aliquot <- colnames(plot_data_mat)
interaction_celltypes <- rownames(plot_data_mat)

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## make color function for heatmap body color
color_blue <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[2]
color_red <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
color_orange <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[5]
color_yellow <- RColorBrewer::brewer.pal(n = 9, name = "YlGnBu")[1]
summary(as.vector(plot_data_mat))
colors_heatmapbody = colorRamp2(c(-1, 
                                  0, 
                                  1), 
                                c(color_blue, "white", color_red))


# make row split ----------------------------------------------------------
therapy_category_vec <- mapvalues(x = interaction_celltypes, from = summary_df$interaction_celltypes, to = as.vector(summary_df$therapy_category))
therapy_category_vec

# make row annotation -----------------------------------------------------
avg_sig

# make column annotation --------------------------------------------------


# make column split -------------------------------------------------------
sampletypes_vec <- mapvalues(x = ids_aliquot, from = idmetadata_df$Aliquot.snRNA.WU, to = as.vector(idmetadata_df$Sample_Type))

# plot hetamp body --------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat,
                             row_names_side = "left",
                             col = colors_heatmapbody,
                             na_col = color_na, border = "black", row_names_max_width = unit(12, "cm"),
                             # show_row_names = T,
                             # row_split = row_split_factor,
                             # row_title_rot = 0, row_title_gp = gpar(fontsize = 15),
                             # # row_labels = factor_cellgroup,
                             # cluster_row_slices = F, 
                             show_row_dend = F, cluster_rows = T, row_split = therapy_category_vec, row_title_side = "right",
                             cluster_columns = T, column_split = sampletypes_vec,
                             # bottom = colanno_obj, show_column_names = T, column_names_gp = gpar(fontsize = 8),
                             # column_title = NA, 
                             # show_column_dend = F, cluster_columns = F, 
                             # column_order = column_order_vec,
                             show_heatmap_legend = T)
p

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "druggable_interactions.", ".png")
png(file2write, width = 1500, height = 2000, res = 150)
print(p)
dev.off()

file2write <- paste0(dir_out, "druggable_interactions.", ".pdf")
pdf(file2write, width = 15, height = 20)
print(p)
dev.off()

