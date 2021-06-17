# Yige Wu @WashU Jul 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
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
## input mutation table
mut_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/mutation/annotate_cptac_sample_by_pbrm1_bap1_mutation/20210412.v1/PBRM1_BAP1_Mutation_Status_By_Case.20210412.v1.tsv")

# get the aliquot IDs for bulk corresponding to the snRNA aliquots --------
metadata_filtered_df <- metadata_bulk_df %>%
  filter(Set.A == "yes") %>%
  filter(Specimen.Label != "CPT0012090003")
metadata_filtered_df$Histologic_Type <- mapvalues(x = metadata_filtered_df$Case.ID, from = case_clinical_df$Case_ID, to = case_clinical_df$Histologic_Type)
metadata_filtered_df <- metadata_filtered_df %>%
  filter(Histologic_Type == "Clear cell renal cell carcinoma")
## add mutation group
metadata_filtered_df$group <- mapvalues(x = metadata_filtered_df$Case.ID, from = mut_df$Case, to = as.vector(mut_df$mutation_category_sim))
metadata_filtered_df$group[as.vector(metadata_filtered_df$group) == as.vector(metadata_filtered_df$Case.ID)] <- "Non-mutants"
table(metadata_filtered_df$group)
## get aliquout ids for the two groups
ids_exp_group1 <- metadata_filtered_df %>%
  filter(Type == "Tumor") %>%
  filter(group %in% c("BAP1 mutated", "Both mutated")) %>%
  filter(!(Case.ID %in% c("C3L-01287")))
ids_exp_group1 <- ids_exp_group1$Specimen.Label ## 16 samples
ids_exp_group2 <- metadata_filtered_df %>%
  filter(Type == "Tumor") %>%
  filter(!(group %in% c("BAP1 mutated", "Both mutated")))
ids_exp_group2 <- ids_exp_group2$Specimen.Label ## 16 samples

# test by wilcox and return values ----------------------------------------
exp_mat <- exp_df[,c(ids_exp_group1, ids_exp_group2)]
gene_index_vec = exp_df$Index
ids_group1 = ids_exp_group1
ids_group2 = ids_exp_group2
registerDoParallel(cores = 3)
start_time <- Sys.time()
test_list<-foreach(g=exp_df$Index) %dopar% {
  exp_raw_vec1 <- unlist(exp_mat[gene_index_vec == g, ids_group1])
  exp_raw_vec2 <- unlist(exp_mat[gene_index_vec == g, ids_group2])
  exp_vec1 <- exp_raw_vec1[!is.na(exp_raw_vec1)]
  exp_vec2 <- exp_raw_vec2[!is.na(exp_raw_vec2)]
  if (length(exp_vec1) >= 5 & length(exp_vec2) >= 5) {
    meddiff_exp <- median(exp_vec1) - median(exp_vec2)
    meddiff_exp
    stat <- wilcox.test(x = exp_vec1, y = exp_vec2)
    p_val <- stat$p.value
    result_list <- list(c(p_val, meddiff_exp, length(exp_vec1), length(exp_vec2)))
  } else {
    result_list <- list(c(NA, NA, length(exp_vec1), length(exp_vec2)))
  }
  return(result_list)
}
end_time <- Sys.time()
end_time - start_time 
## 0.4593129 secs for 100
## 3.820401 secs for 1000 (2 cores), 1.818039 secs for 1000 (3 cores)
## 23.30472 secs for 11355 (3 cores)

## make test result into a data frame
test_df <- data.frame(matrix(data = unlist(test_list), ncol = 4, byrow = T))
colnames(test_df) <- c("p_val", "meddiff_exp", "number_bap1tumors", "number_other_tumors")
test_df$gene_symbol <- exp_df$Index
## adjust p value
test_df$fdr <- p.adjust(p = test_df$p_val, method = "fdr")

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Bulk_Protein_BAP1_vs_Others.Wilcox.", run_id, ".tsv")
write.table(x = test_df, file = file2write, quote = F, sep = "\t", row.names = F)
