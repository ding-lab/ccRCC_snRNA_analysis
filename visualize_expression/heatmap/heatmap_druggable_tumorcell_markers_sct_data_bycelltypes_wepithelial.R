# Yige Wu @WashU Sep 2020

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
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
avgexp_df <- fread(input = "./Resources/Analysis_Results/average_expression/avgexp_sct_data_bycelltypew_epithelial_katmai/20210709.v1/33_aliquot_merged.avgexp.SCT.data.Cell_group_w_epithelialcelltypes.20210709.v1.tsv", data.table = F)
genes_process_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/filter_all_ccRCC_specific_markers_by_druggability/20220715.v1/ccRCC_markers.Druggable.20220715.v1.tsv")

# specify pairs to filter -------------------------------------------------
genes_filter <- genes_process_df$Gene

# format the column names to only aliquot id ------------------------------
## filtr the rows
plot_data_df <- avgexp_df %>%
  rename(gene = V1) %>%
  filter(gene %in% genes_filter)
## filter the columns and make data matrix
plot_data_raw_mat <- as.matrix(plot_data_df[,-1])
plot_data_mat <- t(apply(plot_data_raw_mat, 1, scale))
rownames(plot_data_mat) <- plot_data_df$gene
colnames(plot_data_mat) <- colnames(plot_data_raw_mat)

# get dimension names -----------------------------------------------------
genes_plot <- rownames(plot_data_mat)
genes_plot
celltypes_plot <- colnames(plot_data_mat)

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## make color function for baseline expression
col_baselineexp <- colorRamp2(c(0, 2), c("white", "orange"))
## make color function for heatmap body colors
summary(as.vector(unlist(plot_data_mat)))
color_blue <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[2]
color_red <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
colors_heatmapbody = colorRamp2(c(-1.5, 
                                  0, 
                                  1.5), 
                                c(color_blue, "white", color_red))
# colors_heatmapbody = colorRamp2(seq(from = 0, to = 9, by = 1), 
#                                 c("white", RColorBrewer::brewer.pal(n = 9, name = "YlOrRd")))
## make colors for bulk protein values
color_orange <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[5]
color_purple <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[4]
colors_bulkpro <- colorRamp2(c(-1, 
                               0, 
                               1), 
                             c(color_purple, "white", color_orange))
## make colors for bulk rna values

# make row split for cell types ------------------------------------------
# col_cellgroups_vec <- mapvalues(x = celltypes_plot, from = celltypes2cellgroup_df$colname_celltype, to = as.vector(celltypes2cellgroup_df$Cell_group.shorter))
# col_cellgroups_factor <- factor(x = col_cellgroups_vec, levels = c("Nephron_Epithelium", "Stroma", "Immune"))

# plot heatmap body -----------------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat, 
                             col = colors_heatmapbody, na_col = color_na, border = "black",
                             ## row
                             show_row_names = T, row_names_side = "left", row_names_gp = gpar(fontface = "italic", fontsize = 15),
                             show_row_dend = F, cluster_rows = T,
                             ## column
                             show_column_names = T, column_names_side = "bottom",
                             column_names_gp = gpar(fontsize = 15),
                             show_column_dend = F, 
                             # column_split = col_cellgroups_factor, cluster_column_slices = F, cluster_columns = F,
                             column_title = NULL,
                             show_heatmap_legend = F)
p
## make legend
list_lgd = list(
  Legend(title = "snRNA expression", title_gp = gpar(fontsize = 12),
         col_fun = colors_heatmapbody, 
         legend_width = unit(3, "cm"),
         direction = "horizontal"))

## save heatmap as png
png(filename = paste0(dir_out, "ligandreceptorgenes", ".png"), 
    width = 600, height = 1000, res = 150)
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()
file2write <- paste0(dir_out, "ligandreceptorgenes", ".pdf")
pdf(file2write, width = 5, height = 12)
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()
## save with no legend
file2write <- paste0(dir_out, "ligandreceptorgenes.nolegend", ".pdf")
pdf(file2write, width = 3, height = 2.5)
draw(object = p)
dev.off()

