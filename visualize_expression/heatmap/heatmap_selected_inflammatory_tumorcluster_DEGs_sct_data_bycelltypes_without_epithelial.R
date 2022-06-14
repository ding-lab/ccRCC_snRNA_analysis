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
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the average expression calculated (SCT)
# avgexp_df <- fread(input = "./Resources/Analysis_Results/average_expression/avgexp_sct_data_bycelltypew_epithelial_katmai/20210709.v1/33_aliquot_merged.avgexp.SCT.data.Cell_group_w_epithelialcelltypes.20210709.v1.tsv", data.table = F)
avgexp_df <- fread(input = "./Resources/Analysis_Results/average_expression/avgexp_sct_data_bycelltypew_epithelial_katmai/20210907.v1/35_aliquot_merged.avgexp.SCT.data.Cell_group_w_epithelialcelltypes.20210907.v1.tsv", data.table = F)

# specify pairs to filter -------------------------------------------------
genes_filter <- c("LILRB1", "IL7R", "IL18R1", "IL18RAP")
genes_filter <- c("LILRB1", "IL7R", "IL15RA")
genes_filter <- c("LILRB1", "B2M")
genes_filter <- c("LILRB1", "IL7", "IL7R", "IL15", "IL15RA", "IL2RB")

# format the column names to only aliquot id ------------------------------
## filtr the rows
plot_data_df <- avgexp_df %>%
  rename(gene = V1) %>%
  filter(gene %in% genes_filter)
## filter the columns and make data matrix
plot_data_raw_mat <- as.matrix(plot_data_df[,paste0("SCT.", c("Myofibroblasts", "Endothelial.cells", "DC", "CD4..T.cells", "B.cells", "Macrophages", "NK.cells", "CD8..T.cells", "Tumor.cells", "Fibroblasts"))])
## add row names
rownames(plot_data_raw_mat) <- plot_data_df$gene
## filter rows based on the average expression
genes_plot <- rownames(plot_data_raw_mat)
genes_plot <- genes_filter[genes_filter %in% genes_plot]
plot_data_mat <- t(apply(plot_data_raw_mat, 1, scale))
rownames(plot_data_mat) <- rownames(plot_data_raw_mat); plot_data_mat <- plot_data_mat[genes_plot,]
colnames(plot_data_mat) <- colnames(plot_data_raw_mat)
## filter column
colnames_plot <- colnames(plot_data_mat)
colnames_plot <- colnames_plot[!(grepl(x = colnames_plot, pattern = "Unknown"))]
colnames_plot <- colnames_plot[!(grepl(x = colnames_plot, pattern = "others"))]
plot_data_mat <- plot_data_mat[,colnames_plot]

# process column name labels ----------------------------------------------
cellgroup_label_df <- data.frame(cell_type13 = c("B-cells", "CD4+ T-cells", "CD8+ T-cells", "DC", "Endothelial cells", "Fibroblasts", "Immune others", "Macrophages", "NK cells", 
                                                 "Normal epithelial cells", "Tumor cells", "Unknown",
                                                 "Proximal tubule", "Loop of Henle", "Distal convoluted tubule", 'Principle cells', "Intercalated cells", "Podocytes"))
cellgroup_label_df <- cellgroup_label_df %>%
  mutate(cell_type13.columnname = gsub(x = cell_type13, pattern = "\\-|\\+| ", replacement = "."))
celltypes_plot <- gsub(x = colnames_plot, pattern = "SCT\\.", replacement = "")
celltypelabels_plot <- mapvalues(x = celltypes_plot, from = cellgroup_label_df$cell_type13.columnname, to = as.vector(cellgroup_label_df$cell_type13))

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## make color function for baseline expression
col_baselineexp <- colorRamp2(c(0, 2), c("white", "orange"))
## make color function for heatmap body colors
summary(as.vector(unlist(plot_data_mat)))
color_blue <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[2]
colors_heatmapbody = colorRamp2(c(-1.5, 
                                  0, 
                                  seq(from = 0.5, to = 2, by = 0.5)), 
                                c(color_blue, "white", RColorBrewer::brewer.pal(n = 4, name = "YlOrRd")))
## make colors for the original unscaled expression
# summary(orig_avgexp_vec)
colors_unscaledexp = circlize::colorRamp2(seq(from = 0, to = 4, by = 0.5), 
                                          RColorBrewer::brewer.pal(name = "RdPu", n = 9))
# make row annotation -----------------------------------------------------
## annotate unscaled expression
orig_avgexp_vec <- rowMeans(x = plot_data_raw_mat, na.rm = T)
## get mean SD
# sd_vec <- mapvalues(x = genes_plot, from = sd_bygene_mean_df$V1, to = as.vector(sd_bygene_mean_df$sd_bysample_mean)); sd_vec <- as.numeric(sd_vec)
## make row annotation
row_anno_obj <- rowAnnotation(#ccRCCvsNontumor_snRNA = anno_simple(x = log2fc_tumorvsnontumor_vec, col = colors_snRNA_fc),
  #                               ccRCCvsPT_snRNA = anno_simple(x = log2fc_tumorvspt_vec, col = colors_snRNA_fc), 
  #                               ccRCCvsNAT_bulkRNA = anno_simple(x = log2fc_bulkrna_vec, col = colors_snRNA_fc), 
  #                               ccRCCvsNAT_bulkProtein = anno_simple(x = log2fc_bulkpro_vec, col = colors_snRNA_fc), 
  #                               is_druggable = anno_simple(x = isdruggable_vec, col = colors_yesno), 
  unscaled_snRNA_exp = anno_simple(x = orig_avgexp_vec, col = colors_unscaledexp),
  #                               # SD_tumorcluster_snRNA_exp = anno_simple(x = sd_vec, col = colors_sd),
  annotation_name_side = "bottom", annotation_name_gp = gpar(fontsize = 10), annotation_width = unit(4, "mm"))
# annotation_name_side = "bottom", annotation_name_gp = gpar(fontsize = 10), annotation_width = unit(20, "mm"))


# make row split for cell types ------------------------------------------
# col_cellgroups_vec <- mapvalues(x = celltypes_plot, from = celltypes2cellgroup_df$colname_celltype, to = as.vector(celltypes2cellgroup_df$Cell_group.shorter))
# col_cellgroups_factor <- factor(x = col_cellgroups_vec, levels = c("Nephron_Epithelium", "Stroma", "Immune"))

# plot heatmap body -----------------------------------------------------------------
# p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat, 
#                              col = colors_heatmapbody, na_col = color_na, border = "black",
#                              ## row
#                              show_row_names = T, row_names_side = "left", row_names_gp = gpar(fontface = "italic", fontsize = 10),
#                              show_row_dend = F, cluster_rows = F,
#                              ## column
#                              show_column_names = T, column_names_side = "bottom",
#                              column_names_gp = gpar(fontsize = 10), column_labels = celltypelabels_plot,
#                              show_column_dend = F, 
#                              # column_split = col_cellgroups_factor, cluster_column_slices = F, cluster_columns = F,
#                              column_title = NULL,
#                              show_heatmap_legend = F)
p <- ComplexHeatmap::Heatmap(matrix = t(plot_data_mat), 
                             col = colors_heatmapbody, na_col = color_na, border = "black",
                             show_column_names = T, column_names_side = "top", column_names_gp = gpar(fontface = "italic", fontsize = 10),
                             show_column_dend = F, cluster_columns = F,
                             show_row_names = T, row_names_side = "left",
                             row_names_gp = gpar(fontsize = 10), row_labels = celltypelabels_plot,
                             show_row_dend = F, 
                             # column_split = col_cellgroups_factor, cluster_column_slices = F, cluster_columns = F,
                             column_title = NULL,
                             show_heatmap_legend = F)
p
## make legend
list_lgd = list(
  Legend(title = "Scaled snRNA\nexpression", title_gp = gpar(fontsize = 10),
         col_fun = colors_heatmapbody, 
         legend_width = unit(2, "cm"),
         direction = "vertical"))

## save heatmap as png
png(filename = paste0(dir_out, paste0(genes_plot, collapse = "_"), ".png"), 
    width = 650, height = 1000, res = 150)
draw(object = p, 
     annotation_legend_side = "right", annotation_legend_list = list_lgd)
dev.off()
file2write <- paste0(dir_out, paste0(genes_plot, collapse = "_"), ".pdf")
# pdf(file2write, width = 4, height = 6.5)
pdf(file2write, width = 3.5, height = 3.2)
draw(object = p, 
     annotation_legend_side = "right", annotation_legend_list = list_lgd)
dev.off()
# ## save with no legend
# file2write <- paste0(dir_out, "ligandreceptorgenes.nolegend", ".pdf")
# pdf(file2write, width = 3, height = 2.5)
# draw(object = p)
# dev.off()

