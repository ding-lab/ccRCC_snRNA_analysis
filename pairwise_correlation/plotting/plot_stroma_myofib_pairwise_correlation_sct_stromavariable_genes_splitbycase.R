# Yige Wu @WashU May 2020
## running on local
## for plotting the aliquot-pairwise correlation coefficients for averaged expression of all genes

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
## input id meta data
id_metadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200505.v1/meta_data.20200505.v1.tsv", data.table = F)
## input the spearman pairwise correlation result
pearson_coef.tumorcellvariable_genes.df <- fread(input = "./Resources/Analysis_Results/pairwise_correlation/calculate_myofibroblasts_pairwise_correlation_sct_stromacellvariable_genes/20200901.v1/pearson_coef.stromavaraible_genes.myofibroblasts.avg_exp.sct.byaliquot.20200901.v1.tsv", data.table = F)
## input barcode 2 cell type
barcod2celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_corrected_by_individual_sample_inspection/20200811.v1/31Aliquot.Barcode2CellType.20200811.v1.tsv")
## specify cell type
celltype_show <- "Myofibroblasts"
## make function
number2rangetext = function(x) {
  if (x < 100) {
    range_low <- floor(x/20)*20
    range_high <- ceiling((x+1)/20)*20
    text_range <- paste0("[", range_low, ",", range_high, ")")
  } else {
    text_range <- ">100"
  }
  return(text_range)
}

# make data matrix for heatmap body ---------------------------------------
## reformat data frame to matrix
plot_data_df <- pearson_coef.tumorcellvariable_genes.df
plot_data_mat <- as.matrix(plot_data_df[,-1])
plot_data_mat %>% head()
rownames(plot_data_mat) <- plot_data_df$V1
### get case name
aliquot_celltype <- rownames(plot_data_mat)
aliquot_ids <- str_split_fixed(string = aliquot_celltype, pattern = "_", n = 2)[,1]
case_ids <- mapvalues(x = aliquot_ids, from = id_metadata_df$Aliquot.snRNA, to = as.vector(id_metadata_df$Case))
ids_aliquot_wu <- mapvalues(x = aliquot_ids, from = id_metadata_df$Aliquot.snRNA, to = as.vector(id_metadata_df$Aliquot.snRNA.WU))
## get suffixes of the aliquot ids
suffixes_aliquot_id <- str_split_fixed(string = ids_aliquot_wu, pattern = "-", n = 3)[,3]
suffixes_aliquot_id

# make row and column split -----------------------------------------------
row_split_vec <- case_ids
length(unique(case_ids))
row_split_factor <- factor(x = row_split_vec, levels = c("C3L-00026", "C3N-00733", "C3L-00448", "C3L-00813", "C3L-00610", 
                                                         "C3N-00495", "C3L-00790", "C3L-00004", "C3L-00583", "C3L-00096", 
                                                         "C3L-00088", "C3L-00416", "C3N-01200", "C3L-00908", "C3L-00010", 
                                                         "C3N-00177", "C3L-00079", "C3L-00917", "C3N-00242", "C3L-01313", 
                                                         "C3N-01213", "C3L-01302", "C3L-01287", "C3L-00359"))
column_split_factor <- row_split_factor

# make row order and column order ----------------------------------------------------------
row_order_vec <- rownames(plot_data_mat)[order(suffixes_aliquot_id)]
column_order_vec <- colnames(plot_data_mat)[order(suffixes_aliquot_id)]

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## make color function for heatmap body colors
col_fun = colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red"))
## make colors for the discrete ranges
colors_numbercellrange_vec <- RColorBrewer::brewer.pal(n = 6, name = "PuBuGn")
names(colors_numbercellrange_vec) <- sapply(X = seq(from = 1, to = 101, by = 20), FUN = number2rangetext)

# make row annotation -----------------------------------------------------
barcod2celltype_df$Id_Case <- mapvalues(x = barcod2celltype_df$orig.ident, from = id_metadata_df$Aliquot.snRNA, to = as.vector(id_metadata_df$Case))
barcod2celltype_df$Id_Aliquot_WU <- mapvalues(x = barcod2celltype_df$orig.ident, from = id_metadata_df$Aliquot.snRNA, to = as.vector(id_metadata_df$Aliquot.snRNA.WU))
## count the number of stroma cells
number_cells_byaliquot <- barcod2celltype_df %>%
  filter(Cell_type.shorter == celltype_show) %>%
  select(Id_Aliquot_WU) %>%
  table() %>%
  data.frame() %>%
  rename(Id_Aliquot_WU = ".")

number_cells_byaliquot_long_df <- barcod2celltype_df %>%
  filter(Cell_type.shorter == celltype_show) %>%
  select(Id_Aliquot_WU, Cell_type.detailed) %>%
  table() %>%
  data.frame()
number_cells_byaliquot_long_df
## map discrete ranges
number_cells_byaliquot$text_range <- sapply(X = number_cells_byaliquot$Freq, FUN = number2rangetext)
rownames(number_cells_byaliquot) <- number_cells_byaliquot$Id_Aliquot_WU
number_cells_byaliquot <- number_cells_byaliquot[ids_aliquot_wu,]
## make row annotation
row_anno = rowAnnotation(Sample_Type_Suffix = anno_simple(x = suffixes_aliquot_id,
                                                          col = colors_tumor_segments,
                                                          width = unit(7, "mm")),
                         Number_Cells = anno_simple(x = number_cells_byaliquot$text_range,
                                                    col = colors_numbercellrange_vec,
                                                    width = unit(7, "mm")),
                         # Histologic_Type = anno_simple(x = omics_long_df$Histologic_Type,
                         #                               col = colors_hist_type,
                         #                               width = unit(1.5, "cm")),
                         annotation_name_side = "top", annotation_name_gp = gpar(fontsize = 10))

# make column annotation --------------------------------------------------
col_anno = HeatmapAnnotation(Sample_Type_Suffix = anno_simple(x = suffixes_aliquot_id,
                                                              col = colors_tumor_segments,
                                                              height = unit(7, "mm")),
                             Number_Cells = anno_simple(x = number_cells_byaliquot$text_range,
                                                        col = colors_numbercellrange_vec,
                                                        height = unit(7, "mm")),
                             # Histologic_Type = anno_simple(x = omics_long_df$Histologic_Type,
                             #                               col = colors_hist_type,
                             #                               width = unit(1.5, "cm")),
                             annotation_name_side = "right", annotation_name_gp = gpar(fontsize = 10))

# plot heatmap body with white-yellow-red ------------------------------------------------------
## make heatmap
p <- Heatmap(matrix = plot_data_mat, 
             width = unit(nrow(plot_data_mat), "cm"), height = unit(ncol(plot_data_mat), "cm"),
             column_labels = ids_aliquot_wu,
             row_labels = ids_aliquot_wu,
             right_annotation = row_anno,
             show_row_names = F, show_column_names = F, 
             row_split = row_split_factor, cluster_row_slices = F, row_order = row_order_vec,
             show_row_dend = F, row_title_rot = 0, row_title_side = "right", row_title_gp = gpar(fontsize = 25, fontface = "bold"),
             row_gap = unit(0, "mm"),
             top_annotation = col_anno,
             column_split = column_split_factor, cluster_column_slices = F, column_order = column_order_vec,
             show_column_dend = F, column_title_side = "top", column_title_rot = 90, column_title_gp = gpar(fontsize = 25, fontface = "bold"),
             column_gap = unit(0, "mm"),
             border = "grey50",
             col = col_fun, 
             show_heatmap_legend = F)
## make legend for heattmap body
heatmap_lgd = Legend(col_fun = col_fun, 
                     title = paste0("Pearson's coeffcient\n(variably expressed genes\nwithin stroma cells)"), 
                     direction = "vertical")
## make legend for top annotation
annotation_lgd = list(
  heatmap_lgd,
  Legend(title = paste0("Number of ", celltype_show),
         labels = names(colors_numbercellrange_vec),
         legend_gp = gpar(fill = colors_numbercellrange_vec)))
## save heatmap as png
png(filename = paste0(dir_out, "heatmap", ".png"), 
    width = 3000, height = 2400, res = 150)
draw(object = p, 
     annotation_legend_side = "right", annotation_legend_list = annotation_lgd)
dev.off()

## save heatmap as pdf
pdf(file = paste0(dir_out, "heatmap",".pdf"), 
    width = 20, height = 16)
draw(object = p, 
     annotation_legend_side = "right", annotation_legend_list = annotation_lgd)
dev.off()


