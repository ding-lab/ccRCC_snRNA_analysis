# Yige Wu @WashU March 2020
## running on local
## for calculating the aliquot-pairwise correlation coefficients for averaged expression of all genes

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
id_metadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv", data.table = F)
## input the spearman pairwise correlation result
pearson_coef.tumorcellvariable_genes.df <- fread(input = "./Resources/Analysis_Results/integration/30_aliquot_integration/pairwise_correlation/calculate_tumor_manualsubcluster_pairwise_correlation_tumorcellvariable_genes/20200325.v1/avg_exp_by_tumorsubluster.tumorcellvaraible_genes.pearson_coef20200325.v1.tsv", data.table = F)

# make data matrix for heatmap body ---------------------------------------
## reformat data frame to matrix
plot_data_df <- pearson_coef.tumorcellvariable_genes.df
plot_data_mat <- as.matrix(plot_data_df[,-1])
plot_data_mat %>% head()
## add row names
rownames(plot_data_mat) <- plot_data_df$V1
plot_data_mat %>% head()
### get aliquot ids and case ids
tumorsubcluster_ids <- rownames(plot_data_mat)
aliquot_ids <- str_split_fixed(string = tumorsubcluster_ids, pattern = "_", n = 2)[,1]
case_ids <- mapvalues(x = aliquot_ids, from = id_metadata_df$Aliquot.snRNA, to = as.vector(id_metadata_df$Case))

# make row annotation -------------------------------------------
### create row annotation
# row_anno = rowAnnotation(foo = anno_text(case_ids, 
#                                          location = 0.5, just = "center",
#                                          gp = gpar(fill = uniq_case_colors[case_ids], col = "white", border = "black"),
#                                          width = max_text_width(case_ids)*1.2))
# 
# # make bottom column annotation -------------------------------------------
# bottom_col_anno = HeatmapAnnotation(foo = anno_text(case_ids, 
#                                                     location = 0.5, just = "center",
#                                                     gp = gpar(fill = uniq_case_colors[case_ids], col = "white", border = "black"),
#                                                     width = max_text_width(case_ids)*1.2))

# plot pearson pairwise correlation for variably expressed genes within tumor cells ------------------------------------------------------
## make color function for heatmap body colors
col_fun = colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red"))
## make heatmap
p <- Heatmap(matrix = plot_data_mat,
             # column_km = 2, column_km_repeats = 100,
             # row_km = 2, row_km_repeats = 100,
             col = col_fun, 
             # right_annotation = row_anno,
             # bottom_annotation = bottom_col_anno, 
             show_heatmap_legend = F)
p
## make legend for heattmap body
heatmap_lgd = Legend(col_fun = col_fun, title = "Pearson's coeffcient\n(variably expressed genes\nwithin tumor cells)", direction = "vertical")
## make legend for top annotation
annotation_lgd = list(
  heatmap_lgd)

## save heatmap
png(filename = paste0(dir_out, "avg_exp.tumorcellvariable_genes.pearson_coef.heatmap.", run_id, ".png"), 
    width = 4500, height = 4000, res = 150)
### combine heatmap and heatmap legend
draw(object = p, 
     annotation_legend_side = "right", annotation_legend_list = annotation_lgd)
dev.off()
