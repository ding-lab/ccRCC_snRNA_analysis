# Yige Wu @WashU Mar 2021

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
dirname_deg_files <- c("./Resources/snRNA_Processed_Data/Differentially_Expressed_Genes/PBRM1_vs_BAP1_NonMutants_Tumorcells/20210429.v1/")
## input mutation category
mut_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/mutation/annotate_cptac_sample_by_pbrm1_bap1_mutation/20210310.v1/PBRM1_BAP1_Mutation_Status_By_Case.20210310.v1.tsv")

# input files -------------------------------------------------------------
deg_sup_df <- NULL
paths_deg_file <- list.files(path = dirname_deg_files, full.names = T)
for (path_deg_file_tmp in paths_deg_file) {
  deg_df_tmp <- fread(data.table = F, input = path_deg_file_tmp)
  easyid_tmp <- deg_df_tmp$easyid_tumor[1]
  caseid_tmp <- gsub(pattern = "\\-T[1-9]", replacement = "", x = easyid_tmp)
  mut_category_tmp <- mut_df$mutation_category_sim[mut_df$Case == caseid_tmp]
  deg_df_tmp$group1_mut_category <- mut_category_tmp
  deg_sup_df <- rbind(deg_df_tmp, deg_sup_df)
}

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "PBRM1_vs_BAP1_NonMutants_DEGs.", run_id, ".tsv")
write.table(x = deg_sup_df, file = file2write, quote = F, sep = "\t", row.names = F)
