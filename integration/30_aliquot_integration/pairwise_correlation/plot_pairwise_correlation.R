# Yige Wu @WashU March 2020
## running on local
## for calculating the aliquot-pairwise correlation coefficients for averaged expression of all genes

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/aes.R")
library(ComplexHeatmap)
library(circlize)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input id meta data
id_metadata_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/make_meta_data/20191105.v1/meta_data.20191105.v1.tsv", data.table = F)

# plot spearman pairwise correlation for all genes ------------------------------------------------------
## input the spearman pairwise correlation result
spearman_rhos.all_genes.df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/pairwise_correlation/calculate_pairwise_correlation_all_genes/20200310.v1/avg_exp.all_genes.spearman_rho.20200310.v1.tsv", data.table = F)
## reformat data frame to matrix
spearman_rhos.all_genes.mat <- spearman_rhos.all_genes.df[,-1]
spearman_rhos.all_genes.mat <- as.matrix(spearman_rhos.all_genes.mat)
spearman_rhos.all_genes.mat %>% head()
rownames(spearman_rhos.all_genes.mat) <- spearman_rhos.all_genes.df$V1
## make function for colors
col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
## make heatmap
p <- Heatmap(matrix = spearman_rhos.all_genes.mat,
             col = col_fun,
             name = "Spearman's rho (all genes)")
p
## save heatmap
png(filename = paste0(dir_out, "avg_exp.all_genes.spearman_rho.heatmap.", run_id, ".png"), 
    width = 1300, height = 1000, res = 150)
print(p)
dev.off()

# plot pearson pairwise correlation for all genes ------------------------------------------------------
## input the spearman pairwise correlation result
pearson_coef.all_genes.df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/pairwise_correlation/calculate_pairwise_correlation_all_genes/20200310.v1/avg_exp.all_genes.pearson_coef.20200310.v1.tsv", data.table = F)
## reformat data frame to matrix
pearson_coef.all_genes.mat <- pearson_coef.all_genes.df[,-1]
pearson_coef.all_genes.mat <- as.matrix(pearson_coef.all_genes.mat)
pearson_coef.all_genes.mat %>% head()
rownames(pearson_coef.all_genes.mat) <- pearson_coef.all_genes.df$V1
## make function for colors
col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
## make heatmap
p <- Heatmap(matrix = pearson_coef.all_genes.mat,
             col = col_fun,
             name = "Pearson's coeffcient\n(all genes)")
p
## save heatmap
png(filename = paste0(dir_out, "avg_exp.all_genes.pearson_coef.heatmap.", run_id, ".png"), 
    width = 1300, height = 1000, res = 150)
print(p)
dev.off()

# plot pearson pairwise correlation for variably expressed genes within tumor cells ------------------------------------------------------
## input the spearman pairwise correlation result
pearson_coef.tumorcellvariable_genes.df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/pairwise_correlation/calculate_pairwise_correlation_tumorcellvariable_genes/20200310.v1/avg_exp.tumorcellvaraible_genes.pearson_coef.20200310.v1.tsv", data.table = F)
## reformat data frame to matrix
plot_data_df <- pearson_coef.tumorcellvariable_genes.df
plot_data_mat <- as.matrix(plot_data_df[,-1])
plot_data_mat %>% head()
rownames(plot_data_mat) <- plot_data_df$V1
## make function for colors
col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
## make case name as row annotation
### get case name
case_ids <- mapvalues(x = rownames(plot_data_mat), from = id_metadata_df$Aliquot.snRNA, to = as.vector(id_metadata_df$Case))
### get unique color for each case
uniq_fill_colors <- Polychrome::dark.colors(n = 24)
names(uniq_fill_colors) <- unique(case_ids)
### create row annotation
row_anno = rowAnnotation(foo = anno_text(case_ids, 
                                   location = 0.5, just = "center",
                                   gp = gpar(fill = uniq_fill_colors[case_ids], col = "white", border = "black"),
                                   width = max_text_width(case_ids)*1.2))
### create column annotation
col_anno = HeatmapAnnotation(foo = anno_text(case_ids, 
                                   location = 0.5, just = "center",
                                   gp = gpar(fill = uniq_fill_colors[case_ids], col = "white", border = "black"),
                                   width = max_text_width(case_ids)*1.2))

## make heatmap
p <- Heatmap(matrix = plot_data_mat,
             col = col_fun, 
             right_annotation = row_anno, bottom_annotation = col_anno,
             show_heatmap_legend = F)
p
## make legend
lgd = Legend(col_fun = col_fun, title = "Pearson's coeffcient\n(variably expressed genes within tumor cells)", direction = "horizontal")
## save heatmap
png(filename = paste0(dir_out, "avg_exp.tumorcellvariable_genes.pearson_coef.heatmap.", run_id, ".png"), 
    width = 1400, height = 1500, res = 150)
### combine heatmap and heatmap legend
draw(object = p, heatmap_legend_side = "top", heatmap_legend_list = lgd)
dev.off()
