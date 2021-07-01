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
dirname_deg_files <- c("./Resources/snRNA_Processed_Data/Differentially_Expressed_Genes/PBRM1_BAP1_vs_NonMutants_Tumorcells/")
## input mutation category
mut_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/mutation/annotate_cptac_sample_by_pbrm1_bap1_mutation/20210310.v1/PBRM1_BAP1_Mutation_Status_By_Case.20210310.v1.tsv")

# preprocess --------------------------------------------------------------
paths_deg_file <- list.files(path = dirname_deg_files, full.names = T)
easyids <- str_split_fixed(string = paths_deg_file, pattern = "\\/", n = 7)[,7]
easyids <- str_split_fixed(string = easyids, pattern = "\\.", n = 2)[,1]
cases <- gsub(pattern = "\\-T[0-9]", replacement = "", x = easyids)
easyids_process <- easyids[cases %in% mut_df$Case[mut_df$mutation_category_sim %in% c("BAP1 mutated", "Both mutated")]]
easyids_process <- easyids_process[!(easyids_process %in% "C3L-01287-T1")]
paths_file_process <- paths_deg_file[easyids %in% easyids_process]

# input files -------------------------------------------------------------
deg_sup_df <- NULL
for (path_deg_file_tmp in paths_file_process) {
  deg_df_tmp <- fread(data.table = F, input = path_deg_file_tmp)
  easyid_tmp <- deg_df_tmp$easyid_tumor[1]
  caseid_tmp <- gsub(pattern = "\\-T[1-9]", replacement = "", x = easyid_tmp)
  mut_category_tmp <- mut_df$mutation_category_sim[mut_df$Case == caseid_tmp]
  deg_df_tmp$group1_mut_category <- mut_category_tmp
  deg_sup_df <- rbind(deg_df_tmp, deg_sup_df)
}

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "PBRM1_vs_NonMutants_DEGs.", run_id, ".tsv")
write.table(x = deg_sup_df, file = file2write, quote = F, sep = "\t", row.names = F)
