# Yige Wu @WashU Sep 2020

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
## input the bulk mutation data
mut_matrix_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/mutation/generate_bulk_mutation_table/20191031.v1/snRNA_ccRCC_Mutation_Table.20191031.v1.csv")

# annotate samples based on PBRM1 & BAP1 mutation -------------------------
mut_matrix_df[is.na(mut_matrix_df)] <- ""
mut_matrix_df <- mut_matrix_df %>%
  mutate(mutation_category = ifelse(PBRM1 == "" & BAP1 == "", "Neither mutated",
                                    ifelse(PBRM1 != "" & BAP1 != "", "Both mutated",
                                           ifelse(PBRM1 == "" & BAP1 != "", "BAP1 mutated", "PBRM1 mutated"))))


# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "PBRM1_BAP1_Mutation_Status_By_Case.", run_id, ".tsv")
write.table(x = mut_matrix_df, file = file2write, quote = F, sep = "\t", row.names = F)
