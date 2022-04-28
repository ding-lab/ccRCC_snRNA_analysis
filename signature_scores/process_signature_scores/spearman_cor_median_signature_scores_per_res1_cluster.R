# Yige Wu @WashU May 2020
## plot cell type on integration UMAP

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}

# input dependencies ------------------------------------------------------
## input median(?) signature scores per cluster
results_mat <- readRDS(file = "./Resources/Analysis_Results/signature_scores/process_signature_scores/make_median_signature_scores_per_res1_cluster/20220427.v1/Median_signature_scores_per_res1_cluster.20220427.v1.RDS")
## input the gene set auto-correlation results as a measure of consistency
sigCorr_df <- fread(data.table = F, input = "./Resources/Analysis_Results/signature_scores/run_vision/getSignatureAutocorrelation_30ccRCC_tumorcellreclustered/20220411.v1/ccRCC.30ccRCC.TumorCellsReclustered.Vision.SignatureAutocorrelation.20220411.v1.tsv")

# do correlation -----------------------------------------------------------------
sigCorr_filtered_df <- sigCorr_df %>%
  filter(grepl(pattern = "HALLMARK", x = gene_set)) %>%
  filter(gene_set != "HALLMARK_UV_RESPONSE") %>%
  filter(FDR < 0.05 & C > 0.1) 
genesets_test <- sigCorr_filtered_df$gene_set
# sigScores_test_mat <- results_mat[genesets_test,]
sigScores_mat <- results_mat[genesets_test,]
sigScores_scaled_mat <- t(apply(sigScores_mat, 1, scale))
colnames(sigScores_scaled_mat) <- colnames(sigScores_mat)
result <- cor(sigScores_scaled_mat, method = "spearman")
## pearson results are not very telling
# result <- cor(sigScores_test_mat, method = "pearson")

library(corrplot)
corrplot(result, type = "upper", order = "hclust", 
         tl.col = "black", tl.srt = 45)

col<- colorRampPalette(c("blue", "white", "red"))(20)
heatmap(x = result, col = col, symm = TRUE)


summary(as.vector(result))
col_fun = colorRamp2(c(0, 0.5, 1), c("blue", "white", "red"))
Heatmap(matrix = result, col = col_fun)

# write output ------------------------------------------------------------
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
## write output
file2write <- paste0(dir_out, "Median_signature_scores_per_res1_cluster.", run_id, ".RDS")
saveRDS(object = plotdata_wide_df, file = file2write, compress = T)


