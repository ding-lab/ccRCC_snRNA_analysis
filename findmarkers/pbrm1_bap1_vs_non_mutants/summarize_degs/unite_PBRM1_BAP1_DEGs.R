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
dir_deg_file_level1 <- "./Resources/snRNA_Processed_Data/Differentially_Expressed_Genes/"
dirs_deg_files <- c("./Resources/snRNA_Processed_Data/Differentially_Expressed_Genes/Tumorcells_vs_PTcells/",
                    "./Resources/snRNA_Processed_Data/Differentially_Expressed_Genes/PBRM1_BAP1_vs_NonMutants_Tumorcells/",
                    "./Resources/snRNA_Processed_Data/Differentially_Expressed_Genes/BAP1_vs_PBRM1_Mutants_Tumorcells/",
                    "./Resources/snRNA_Processed_Data/Differentially_Expressed_Genes/PBRM1_vs_BAP1_Mutants_Tumorcells/")
## input mutation category
mut_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/mutation/annotate_cptac_sample_by_pbrm1_bap1_mutation/20210310.v1/PBRM1_BAP1_Mutation_Status_By_Case.20210310.v1.tsv")

# input files -------------------------------------------------------------
dirnames_deg_files <- list.dirs(path = dir_deg_file_level1, full.names = F)
dirnames_deg_files <- dirnames_deg_files[grepl(pattern = "vs", x = dirnames_deg_files)]
deg_sup_df <- NULL
for (dirname_deg_file_tmp in dirnames_deg_files) {
  dir_deg_file_tmp <- paste0(dir_deg_file_level1, dirname_deg_file_tmp, "/")
  paths_deg_file <- list.files(path = dir_deg_file_tmp, full.names = T)
  for (path_deg_file_tmp in paths_deg_file) {
    deg_df_tmp <- fread(data.table = F, input = path_deg_file_tmp)
    if (dirname_deg_file_tmp %in% c("BAP1_vs_PBRM1_Mutants_Tumorcells", "PBRM1_vs_BAP1_Mutants_Tumorcells")) {
      deg_df_tmp$comparison <- dirname_deg_file_tmp
    } else {
      easyid_tmp <- deg_df_tmp$easyid_tumor[1]
      caseid_tmp <- gsub(pattern = "\\-T[1-9]", replacement = "", x = easyid_tmp)
      mut_category_tmp <- mut_df$mutation_category_sim[mut_df$Case == caseid_tmp]
      if (dirname_deg_file_tmp == "PBRM1_BAP1_vs_NonMutants_Tumorcells") {
        deg_df_tmp$comparison <- ifelse(mut_category_tmp == "Both mutated", "PBRM1_BAP1_vs_NonMutants_Tumorcells",
                                        ifelse(mut_category_tmp == "BAP1 mutated", "BAP1_vs_NonMutants_Tumorcells", "PBRM1_vs_NonMutants_Tumorcells"))
      }
      if (dirname_deg_file_tmp == "Tumorcells_vs_PTcells") {
        deg_df_tmp$comparison <- ifelse(mut_category_tmp == "Non-mutants", "NonMutant_Tumorcells_vs_PTcells",
                                        ifelse(mut_category_tmp == "Both mutated", "PBRM1_BAP1_Mutant_Tumorcells_vs_PTcells",
                                               ifelse(mut_category_tmp == "PBRM1 mutated", "PBRM1_Tumorcells_vs_PTcells", "BAP1_Tumorcells_vs_PTcells")))
      }
    }
    deg_sup_df <- rbind(deg_df_tmp, deg_sup_df)
  }
}

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "PBRM1_BAP1_DEGs.", run_id, ".tsv")
write.table(x = deg_sup_df, file = file2write, quote = F, sep = "\t", row.names = F)
