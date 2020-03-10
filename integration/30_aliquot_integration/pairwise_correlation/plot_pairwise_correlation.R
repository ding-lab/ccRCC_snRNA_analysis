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

# plot spearman pairwise correlation for all genes ------------------------------------------------------
## input the spearman pairwise correlation result
spearman_rhos.all_genes.df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/pairwise_correlation/calculate_pairwise_correlation_all_genes/20200310.v1/avg_exp.all_genes.spearman_rho.20200310.v1.tsv", data.table = F)
## reformat data frame to matrix
spearman_rhos.all_genes.mat <- spearman_rhos.all_genes.df[,-1]
spearman_rhos.all_genes.mat <- as.matrix(spearman_rhos.all_genes.mat)
spearman_rhos.all_genes.mat %>% head()
rownames(spearman_rhos.all_genes.mat) <- spearman_rhos.all_genes.df$V1
## make function for colors
col_fun = colorRamp2(c(0, 0.5, 1), c("green", "white", "red"))
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
col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "yellow", "red"))
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
