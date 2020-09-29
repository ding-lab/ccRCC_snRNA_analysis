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
summary_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/filter_interactions/filter_stroma_tumor_interactions_bynumbersigcases_topavgsigmean/20200923.v1/cellphonedb.summary_across_samples_by_pair.stroma&tumor.filtered.20200923.v1.tsv")
## input cell-cell interaction values by case
cellphone_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/annotate_interactions/annotate_cellphone_out/20200923.v1/cell.phone.res.total.run20200818.filtered.formatted.txt")

# specify pairs to filter -------------------------------------------------
nrow(summary_df)
pair_cell.types_process <- summary_df$pair_cell.types

# make data matrix for the heatmap body -------------------------------------------
plot_data_df <- cellphone_df %>%
  filter(pair_cell.types %in% pair_cell.types_process)
plot_data_wide_df <- dcast(data = plot_data_df, formula = pair_cell.types ~ Easy_id, value.var = "value")
plot_data_mat <- as.matrix(plot_data_wide_df[, -1])
rownames(plot_data_mat) <- plot_data_wide_df$pair_cell.types

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## make color function for heatmap body colors
color_red <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
color_orange <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[5]
color_yellow <- RColorBrewer::brewer.pal(n = 9, name = "YlGnBu")[1]

# summary(as.vector(plot_data_mat))
colors_heatmapbody = colorRamp2(c(0, 
                                  10, 
                                  20), 
                                c("white", color_orange, color_red))


# plot hetamp body --------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat,
                             row_names_side = "left",
                             col = colors_heatmapbody,
                             na_col = color_na, border = "black", row_names_max_width = unit(10, "cm"),
                             # show_row_names = T,
                             # row_split = row_split_factor,
                             # row_title_rot = 0, row_title_gp = gpar(fontsize = 15),
                             # # row_labels = factor_cellgroup,
                             # cluster_row_slices = F, 
                             show_row_dend = F,
                             # bottom = colanno_obj, show_column_names = T, column_names_gp = gpar(fontsize = 8),
                             # column_title = NA, 
                             # show_column_dend = F, cluster_columns = F, 
                             # column_order = column_order_vec,
                             show_heatmap_legend = T)
p

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "stroma_tumor_interactions.", ".png")
png(file2write, width = 1500, height = 2000, res = 150)
print(p)
dev.off()

file2write <- paste0(dir_out, "stroma_tumor_interactions.", ".pdf")
pdf(file2write, width = 15, height = 20)
print(p)
dev.off()

