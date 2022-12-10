# Yige Wu @WashU June 2022

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "data.table",
  "stringr"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    install.packages(pkg_name_tmp, dependencies = T)
  }
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg_name_tmp)
  }
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input data --------------------------------------------------------------
## input bulk gene expression
exp_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/mRNA/CPTAC_ccRCC/CPTAC_ccRCC_discovery_tumor_mRNA_FPKM_UQ_log2_v1.0.tsv")
## input the tumor-cell-intrinsic inflammation signature
markers_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/filter_degs/overlap_intrapatienttumorclusters_selected_top_vs_bottom_with_tumor_markers/20221207.v1/selectedpathway_markers_overlap_tumorcell_markers.20221207.v1.tsv")
## input bulk meta data
metadata_bulk_df <- fread("./Resources/Bulk_Processed_Data/Case_ID/CPTAC_ccRCC_discovery_caseID_v1.0.tsv")

# preprocess ------------------------------------------------------
metadata_filtered_df <- metadata_bulk_df %>%
  filter(Histologic_Type == "Clear cell renal cell carcinoma")


# process by gene set -----------------------------------------------------
scores_all_df <- NULL
for (geneset_tmp in unique(markers_df$geneset)) {
  exp_markers_df <- exp_df %>%
    filter(gene_name %in% markers_df$genesymbol_deg[markers_df$geneset == geneset_tmp])
  
  exp_markers_mat <- exp_markers_df[,3:ncol(exp_markers_df)]
  colnames(exp_markers_mat) <- gsub(x = colnames(exp_markers_mat), pattern = "\\-T", replacement = "")
  exp_markers_mat <- exp_markers_mat[, metadata_filtered_df$CASE_ID]
  exp_markers_scale_mat <- t(apply(exp_markers_mat, 1, scale))
  rownames(exp_markers_scale_mat) <- exp_markers_df$gene_name
  colnames(exp_markers_scale_mat) <- colnames(exp_markers_mat)
  scores_vec <- colMeans(exp_markers_scale_mat, na.rm = T)*100
  scores_df <- data.frame(case = names(scores_vec), geneset_score = scores_vec, geneset = geneset_tmp)
  scores_all_df <- rbind(scores_all_df, scores_df)
}

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "selectedpathway_tumorcell_signature_scores.", run_id, ".tsv")
write.table(x = scores_all_df, file = file2write, quote = F, sep = "\t", row.names = F)

