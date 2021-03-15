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
## input the bulk mutation data
maf_df <- fread(data.table = F, input = "~/Box/Ding_Lab/Databases/cBioPortal/CCRCC/kirc_tcga_pub/data_mutations_extended.txt")
metadata_df <- fread(data.table = F, input = "~/Box/Ding_Lab/Databases/cBioPortal/CCRCC/kirc_tcga_pub/data_clinical_sample.txt")
metadata_df <- fread(data.table = F, input = "../ccRCC_snATAC/Resources/Knowledge/TCGA_ATAC/TCGA_ATAC_peak.all.probeMap")

# annotate samples based on PBRM1 & BAP1 mutation -------------------------
mut_matrix_wide_df <- get_mutation_aachange_matrix(maf = maf_df, pair_tab = ccRCC_SMGs)
mut_matrix_df <- as.data.frame(t(mut_matrix_wide_df[,-1]))
mut_matrix_df$Case <- rownames(mut_matrix_df)
mut_matrix_df[is.na(mut_matrix_df)] <- ""
mut_matrix_df <- mut_matrix_df %>%
  mutate(mutation_category_sim = ifelse(PBRM1 == "" & BAP1 == "", "Non-mutants",
                                        ifelse(PBRM1 != "" & BAP1 != "", "Both mutated",
                                               ifelse(PBRM1 == "" & BAP1 != "", "BAP1 mutated", "PBRM1 mutated")))) %>%
  mutate(mutation_category = ifelse(PBRM1 == "" & BAP1 == "",
                                    ifelse(KDM5C == "" & SETD2 == "", "Non-mutants", "Others"),
                                    ifelse(PBRM1 != "" & BAP1 != "", "Both mutated",
                                           ifelse(PBRM1 == "" & BAP1 != "", "BAP1 mutated", "PBRM1 mutated"))))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "PBRM1_BAP1_Mutation_Status_By_Case.", run_id, ".tsv")
write.table(x = mut_matrix_df, file = file2write, quote = F, sep = "\t", row.names = F)
