# Yige Wu @WashU March 2020
## running on local
## for calculating the correlation of average expression between tumor subclusters

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the variable gene list
variable_genes_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/findvariablefeatures/findvariablefeatures_tumor_cells/20200310.v1/findvariablefeatures_tumor_cells.20200310.v1.tsv", data.table = F)
## input the average expression calculated (RNA)
avg.exp.mat <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/averageexpression/averageexpression_tumor_cells_by_manual_subcluster/20200325.v1/averageexpression_tumor_cells_by_manual_subcluster.20200325.v1.tsv", data.table = F)

# format the column names to only aliquot id ------------------------------
avg.exp.mat <- avg.exp.mat %>%
  rename(gene = V1)
data_col_names <- colnames(avg.exp.mat)[-1]
data_col_names.changed <- str_split_fixed(string = data_col_names, pattern = "\\.", n = 2)[,2]
## rename the data frame
colnames(avg.exp.mat) <- c("gene", data_col_names.changed)
## remove the subcluster without NA as manual subcluster id
data_col_names.filtered <- data_col_names.changed[!grepl(pattern = "MCNA", x = data_col_names.changed)]
avg.exp.mat <- avg.exp.mat[, c("gene", data_col_names.filtered)]
## take out the genes with 0 average expression
avg.exp.mat <- avg.exp.mat[rowSums(avg.exp.mat[,-1]) > 0, ]

# run and write spearman pairwise correlation ------------------------------------------------
## For cor(), if method is "kendall" or "spearman", Kendall's tau or Spearman's rho statistic is used to estimate a rank-based measure of association. These are more robust and have been recommended if the data do not necessarily come from a bivariate normal distribution.
spearman_rhos <- cor(x = avg.exp.mat[, -1], use="complete.obs", method = "spearman")
spearman_rhos
write.table(x = spearman_rhos, file = paste0(dir_out, "avg_exp_by_tumorsubluster.", "all_genes.", "spearman_rho.", run_id, ".tsv"), quote = F, sep = "\t", row.names = T)

# run and write pearson pairwise correlation ------------------------------------------------
pearson_coef <- cor(x = avg.exp.mat[, -1], use="complete.obs", method = "pearson")
pearson_coef
write.table(x = pearson_coef, file = paste0(dir_out, "avg_exp_by_tumorsubluster.", "all_genes.", "pearson_coef", run_id, ".tsv"), quote = F, sep = "\t", row.names = T)