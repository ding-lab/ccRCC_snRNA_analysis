# Yige Wu @WashU Apr 2022
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
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
results_df <- fread(input = "./Resources/Analysis_Results/signature_scores/compare_signature_scores/wilcox_pairwise_bynewcluster_bygeneset_30ccRCC_tumorcellsreclustered/20220407.v1/wilcox.pairwise.bynewcluster.30ccRCCtumorcellreclustered,20220407.v1.tsv", data.table = F)

# process -----------------------------------------------------------------
results_filtered_df <- results_df %>%
  filter(ident.1 == "C2" & ident.2 == "C1")
results_filtered_df$fdr <- p.adjust(p = results_filtered_df$p_val, method = "fdr")
results_filtered_df2 <- results_filtered_df %>%
  filter(fdr < 0.01)
