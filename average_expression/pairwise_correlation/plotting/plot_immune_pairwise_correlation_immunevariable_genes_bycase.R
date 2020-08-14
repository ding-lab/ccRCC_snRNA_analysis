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
pearson_coef.tumorcellvariable_genes.df <- fread(input = "./Resources/Analysis_Results/average_expression/pairwise_correlation/calculate_immune_pairwise_correlation_immunecellvariable_genes/20200811.v1/avg_exp_bycellgroup_byaliquot.immune.immunevaraible_genes.pearson_coef20200811.v1.tsv", data.table = F)
## input barcode 2 cell type
barcod2celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_corrected_by_individual_sample_inspection/20200811.v1/31Aliquot.Barcode2CellType.20200811.v1.tsv")
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
aliquot_wu_ids <- mapvalues(x = aliquot_ids, from = id_metadata_df$Aliquot.snRNA, to = as.vector(id_metadata_df$Aliquot.snRNA.WU))

# make row annotation -----------------------------------------------------
barcod2celltype_df$Id_Case <- mapvalues(x = barcod2celltype_df$orig.ident, from = id_metadata_df$Aliquot.snRNA, to = as.vector(id_metadata_df$Case))
barcod2celltype_df$Id_Aliquot_WU <- mapvalues(x = barcod2celltype_df$orig.ident, from = id_metadata_df$Aliquot.snRNA, to = as.vector(id_metadata_df$Aliquot.snRNA.WU))
number_cells_byaliquot <- barcod2celltype_df %>%
  filter(Cell_group.shorter == "Immune") %>%
  select(Id_Aliquot_WU) %>%
  table() %>%
  data.frame() %>%
  rename(Id_Aliquot_WU = ".")
## map discrete ranges
number_cells_byaliquot$text_range <- sapply(X = number_cells_byaliquot$Freq, FUN = number2rangetext)
rownames(number_cells_byaliquot) <- number_cells_byaliquot$Id_Aliquot_WU
number_cells_byaliquot <- number_cells_byaliquot[aliquot_wu_ids,]
## make colors for the discrete ranges
colors_numbercellrange_vec <- RColorBrewer::brewer.pal(n = 6, name = "BuPu")
names(colors_numbercellrange_vec) <- sapply(X = seq(from = 1, to = 101, by = 20), FUN = number2rangetext)
## make row annotation
row_left_anno = rowAnnotation(Number_StromaCells = anno_simple(x = number_cells_byaliquot$text_range,
                                                               col = colors_numbercellrange_vec,
                                                               width = unit(5, "mm")))

# plot heatmap body with white-yellow-red ------------------------------------------------------
## make color function for heatmap body colors
col_fun = colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red"))
## make heatmap
p <- Heatmap(matrix = plot_data_mat, 
             column_labels = aliquot_wu_ids,
             row_labels = aliquot_wu_ids,
             right_annotation = row_left_anno,
             show_row_names = T, show_column_names = T, 
             row_split = case_ids, cluster_row_slices = F,
             show_row_dend = F, row_title_rot = 0, row_title_gp = gpar(fontsize = 10, fontface = "bold"),
             row_gap = unit(0, "mm"),
             column_split = case_ids, cluster_column_slices = F, column_title_side = "top",
             show_column_dend = F, column_title_rot = 90, column_title_gp = gpar(fontsize = 10, fontface = "bold"),
             column_gap = unit(0, "mm"),
             border = "grey50",
             col = col_fun, 
             show_heatmap_legend = F)
p
## make legend for heattmap body
heatmap_lgd = Legend(col_fun = col_fun, 
                     title = "Pearson's coeffcient\n(variably expressed genes\nwithin immune cells)", 
                     direction = "vertical")
## make legend for top annotation
annotation_lgd = list(
  heatmap_lgd,
  Legend(title = "Number of Immune Cells",
         labels = names(colors_numbercellrange_vec),
         legend_gp = gpar(fill = colors_numbercellrange_vec)))
## save heatmap as png
png(filename = paste0(dir_out, "avg_exp.tumorcellvariable_genes.pearson_coef.heatmap.WYR.", run_id, ".png"), 
    width = 1500, height = 1200, res = 150)
draw(object = p, 
     annotation_legend_side = "right", annotation_legend_list = annotation_lgd)
dev.off()

## save heatmap as pdf
pdf(file = paste0(dir_out, "avg_exp.tumorcellvariable_genes.pearson_coef.heatmap.WYR.", run_id, ".pdf"), 
    width = 8, height = 6)
draw(object = p, 
     annotation_legend_side = "right", annotation_legend_list = annotation_lgd)
dev.off()


