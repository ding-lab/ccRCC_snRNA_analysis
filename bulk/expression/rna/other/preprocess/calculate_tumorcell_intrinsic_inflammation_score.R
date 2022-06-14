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
inflamm_markers_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/filter_degs/overlap_intrapatienttumorclusters_inflammatory_top_vs_bottom_with_tumor_markers/20220612.v1/inflammatory_tumorcell_markers.20220612.v1.tsv")
## input bulk meta data
metadata_bulk_df <- fread("./Resources/Bulk_Processed_Data/Case_ID/CPTAC_ccRCC_discovery_caseID_v1.0.tsv")

# preprocess ------------------------------------------------------
metadata_filtered_df <- metadata_bulk_df %>%
  filter(Histologic_Type == "Clear cell renal cell carcinoma")
exp_inflamm_markers_df <- exp_df %>%
  filter(gene_name %in% inflamm_markers_df$genesymbol_deg)

# process -----------------------------------------------------------------
exp_inflamm_markers_mat <- exp_inflamm_markers_df[,3:ncol(exp_inflamm_markers_df)]
colnames(exp_inflamm_markers_mat) <- gsub(x = colnames(exp_inflamm_markers_mat), pattern = "\\-T", replacement = "")
exp_inflamm_markers_mat <- exp_inflamm_markers_mat[, metadata_filtered_df$CASE_ID]
exp_inflamm_markers_scale_mat <- t(apply(exp_inflamm_markers_mat, 1, scale))
rownames(exp_inflamm_markers_scale_mat) <- exp_inflamm_markers_df$gene_name
colnames(exp_inflamm_markers_scale_mat) <- colnames(exp_inflamm_markers_mat)
inflamm_scores_vec <- colMeans(exp_inflamm_markers_scale_mat, na.rm = T)*100
summary(inflamm_scores_vec)
names(inflamm_scores_vec)
inflamm_scores_df <- data.frame(case = names(inflamm_scores_vec), tumorcellintrinsic_inflamm_score = inflamm_scores_vec)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "tumorcell_intrinsic_inflammation_signature_scores.", run_id, ".tsv")
write.table(x = inflamm_scores_df, file = file2write, quote = F, sep = "\t", row.names = F)

