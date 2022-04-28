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
results_df <- fread(data.table = F, input = "./Resources/Analysis_Results/signature_scores/compare_signature_scores/wilcox_eachcluster_vs_rest_bygeneset_res1_30ccRCC_tumorcellsreclustered/20220411.v2/wilcox.eachcluster_vs_rest.res1.30ccRCCtumorcellreclustered,20220411.v2.tsv")

# make matrix data for heatmap body color-----------------------------------------------------------------
## extract the data for the matrix
plotdata_df <- results_df %>%
  mutate(x_plot = gene_set) %>%
  mutate(y_plot = cluster) %>%
  mutate(value = median_sigScore)
plotdata_wide_df <- dcast(data = plotdata_df, formula = x_plot ~ y_plot, value.var = "value")
colnames(plotdata_wide_df)[1] <- "gene_set"
plotdata_mat <- as.matrix(plotdata_wide_df[,-1])
rownames(plotdata_mat) <- plotdata_wide_df$gene_set

# write output ------------------------------------------------------------
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
## write output
file2write <- paste0(dir_out, "Median_signature_scores_per_res1_cluster.", run_id, ".RDS")
saveRDS(object = plotdata_mat, file = file2write, compress = T)


