# Yige Wu @WashU Jul 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
exp_df <- fread("./Resources/Bulk_Processed_Data/Protein/6_CPTAC3_CCRCC_Whole_abundance_gene_protNorm=2_CB.tsv", data.table = F)
## input bulk meta data
metadata_bulk_df <- fread("./Resources/Bulk_Processed_Data/Meta_Data/cptac-metadata.csv")
## input clinical info
case_clinical_df <- readxl::read_excel(path = "./Resources/Bulk_Processed_Data/CPTAC3-ccRCC-SupplementaryTables_Final/Table S1.xlsx", sheet = "ccrcc_clinical_characteristics")
  
# get the aliquot IDs for bulk corresponding to the snRNA aliquots --------
metadata_filtered_df <- metadata_bulk_df %>%
  filter(Set.A == "yes") %>%
  filter(Specimen.Label != "CPT0012090003")
metadata_filtered_df$Histologic_Type <- mapvalues(x = metadata_filtered_df$Case.ID, from = case_clinical_df$Case_ID, to = case_clinical_df$Histologic_Type)
metadata_filtered_df <- metadata_filtered_df %>%
  filter(Histologic_Type == "Clear cell renal cell carcinoma")
ids_exp_normal <- metadata_filtered_df$Specimen.Label[metadata_filtered_df$Type == "Normal"]
ids_exp_normal

case_ids2plot <- mapvalues(x = ids_exp_normal, from = metadata_filtered_df$Specimen.Label, to = metadata_filtered_df$Case.ID)
case_ids2plot

ids_exp_tumor <- mapvalues(x = case_ids2plot, from = metadata_filtered_df$Case.ID[metadata_filtered_df$Type == "Tumor"], to = as.vector(metadata_filtered_df$Specimen.Label[metadata_filtered_df$Type == "Tumor"]))
ids_exp_tumor

# test by wilcox and return values ----------------------------------------
test_list <- lapply(exp_df$Index, function(g, exp_mat, gene_index_vec) {
  exp_raw_t_vec <- unlist(exp_mat[gene_index_vec == g, ids_exp_tumor])
  exp_raw_n_vec <- unlist(exp_mat[gene_index_vec == g, ids_exp_normal])
  index_nonna <- (!is.na(exp_raw_t_vec) & !is.na(exp_raw_n_vec))
  exp_t_vec <- exp_raw_t_vec[index_nonna]
  exp_n_vec <- exp_raw_n_vec[index_nonna]
  if (length(exp_t_vec) >= 5) {
    meddiff_exp <- median(exp_t_vec) - median(exp_n_vec)
    meddiff_exp
    stat <- wilcox.test(x = exp_t_vec, y = exp_n_vec, paired = T)
    p_val <- stat$p.value
    result_list <- list(c(p_val, meddiff_exp, length(exp_t_vec)))
  } else {
    result_list <- list(c(NA, NA, length(exp_t_vec)))
  }
  return(result_list)
}, exp_mat = exp_df[,c(ids_exp_tumor, ids_exp_normal)], gene_index_vec = exp_df$Index)
## make test result into a data frame
test_df <- data.frame(matrix(data = unlist(test_list), ncol = 3, byrow = T))
colnames(test_df) <- c("p_val", "meddiff_exp", "number_values")
test_df$gene_symbol <- exp_df$Index
## adjust p value
test_df$fdr <- p.adjust(p = test_df$p_val, method = "fdr")

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Bulk_Protein_Tumor_vs_Normal.Wilcox.", run_id, ".tsv")
write.table(x = test_df, file = file2write, quote = F, sep = "\t", row.names = F)
