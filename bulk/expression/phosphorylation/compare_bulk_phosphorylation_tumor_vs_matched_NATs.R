# Yige Wu @WashU Jul 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(doParallel)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
exp_df <- fread(data.table = F, input = "~/Box/CPTAC_ccRCC/Data_Freeze_1.1/CPTAC_ccRCC_Discovery/Phosphoproteome/TMT/6_CPTAC3_CCRCC_Phospho_abundance_phosphosite_protNorm=2_CB_imputed.tsv")
metadata_df <- fread(data.table = F, input = "~/Box/CPTAC_ccRCC/Data_Freeze_1.1/CPTAC_ccRCC_Discovery/Case_ID/CPTAC_ccRCC_discovery_caseID_v1.0.tsv")

# get the aliquot IDs for bulk corresponding to the snRNA aliquots --------
metadata_filtered_df <- metadata_df %>%
  filter(Specimen.Label.tumor != "CPT0012090003") %>%
  filter(Histologic_Type == "Clear cell renal cell carcinoma") %>%
  filter(!is.na(Specimen.Label.normal))
## get aliquout ids for the two groups
ids_exp_group1 <- metadata_filtered_df$Specimen.Label.tumor; length(ids_exp_group1) ## 81 samples
ids_exp_group2 <- metadata_filtered_df$Specimen.Label.normal
ids_exp_group2 <- ids_exp_group2[ids_exp_group2 %in% colnames(exp_df)]; length(ids_exp_group2) ## 81 samples
## process phosphorylation data
phosphosite_idx_df <- data.frame(SUBSTRATE = exp_df$Gene, id_phosphosite = exp_df$Index)
phosphosite_idx_df <- phosphosite_idx_df %>%
  mutate(SUB_MOD_RSD = str_split_fixed(string = id_phosphosite, pattern = "_", n = 7)[,7]) %>%
  filter(SUB_MOD_RSD != "")

# test by wilcox and return values ----------------------------------------
exp_mat <- exp_df[exp_df$Index %in% phosphosite_idx_df$id_phosphosite,c(ids_exp_group1, ids_exp_group2)]
gene_index_vec = phosphosite_idx_df$id_phosphosite
ids_group1 = ids_exp_group1
ids_group2 = ids_exp_group2
registerDoParallel(cores = 3)
start_time <- Sys.time()
test_list<-foreach(g=gene_index_vec) %dopar% {
  exp_raw_vec1 <- unlist(exp_mat[gene_index_vec == g, ids_group1])
  exp_raw_vec2 <- unlist(exp_mat[gene_index_vec == g, ids_group2])
  exp_vec1 <- exp_raw_vec1[!is.na(exp_raw_vec1)]
  exp_vec2 <- exp_raw_vec2[!is.na(exp_raw_vec2)]
  meddiff_exp <- median(exp_vec1) - median(exp_vec2)
  meddiff_exp
  stat <- wilcox.test(x = exp_vec1, y = exp_vec2, paired = T)
  p_val <- stat$p.value
  result_list <- c(p_val, meddiff_exp)
  return(result_list)
}
end_time <- Sys.time()
end_time - start_time 

## make test result into a data frame
test_df <- do.call(rbind.data.frame, test_list)
# test_df <- data.frame(matrix(data = unlist(test_list), ncol = 4, byrow = T))
colnames(test_df) <- c("p_val", "meddiff_exp")
test_df <- cbind(phosphosite_idx_df, test_df)
## adjust p value
test_df$fdr <- p.adjust(p = test_df$p_val, method = "fdr")

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Bulk_Phosphosite_Tumor_vs_NAT.Matched.Wilcox.", run_id, ".tsv")
write.table(x = test_df, file = file2write, quote = F, sep = "\t", row.names = F)

# summarize ---------------------------------------------------------------
test_df %>%
  filter(fdr < 0.05) %>%
  filter(meddiff_exp > 0) %>%
  nrow()

test_df %>%
  nrow()
