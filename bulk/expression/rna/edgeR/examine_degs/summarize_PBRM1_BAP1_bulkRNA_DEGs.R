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
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/edgeR/run_deg_analysis/run_PBRM1_BAP1_deg_on_cptac_ccRCC_discovery_cases/20210329.v1/PBRM1_BAP1_DEGs.glmQLFTest.OutputTables.tsv")

# filter genes ---------------------------------------------
deg_filtered_df <- deg_df %>%
  filter(FDR < 0.05) %>%
  mutate(genesymbol_deg = ifelse(!is.na(hgnc_symbol) & hgnc_symbol != "", hgnc_symbol, clone_based_ensembl_gene)) %>%
  filter(!is.na(genesymbol_deg) & genesymbol_deg != "") %>%
  mutate(row_id = paste0(comparison, "-", genesymbol_deg)) %>%
  filter(!duplicated(row_id))
deg_wide_df <- dcast(data = deg_filtered_df, formula = genesymbol_deg~comparison, value.var = "logFC")

# filter BAP1-specific DEGs ----------------------------------------------
deg_wide_df <- deg_wide_df %>%
  mutate(BAP1_Mutated_vsOthers.bulkRNA = ifelse(BAP1_Mutated_vs_PBRM1_Mutated > 0 & BAP1_Mutated_vs_NonMutants > 0 & BAP1_Mutated_vs_NAT > 0, "Up",
                                        ifelse(BAP1_Mutated_vs_PBRM1_Mutated < 0 & BAP1_Mutated_vs_NonMutants < 0 & BAP1_Mutated_vs_NAT < 0, "Down", "Inconsistent")))
deg_bap1_df <- deg_wide_df %>%
  filter(BAP1_Mutated_vsOthers.bulkRNA != "Inconsistent")
nrow(deg_bap1_df)
table(deg_bap1_df$BAP1_Mutated_vsOthers.bulkRNA)

# filter PBRM1-specific DEGs ----------------------------------------------
deg_wide_df <- deg_wide_df %>%
  mutate(PBRM1_Mutated_vsOthers.bulkRNA = ifelse(BAP1_Mutated_vs_PBRM1_Mutated < 0 & PBRM1_Mutated_vs_NonMutants > 0 & PBRM1_Mutated_vs_NAT > 0, "Up",
                                         ifelse(BAP1_Mutated_vs_PBRM1_Mutated > 0 & PBRM1_Mutated_vs_NonMutants < 0 & PBRM1_Mutated_vs_NAT < 0, "Down", "Inconsistent")))
deg_pbrm1_df <- deg_wide_df %>%
  filter(PBRM1_Mutated_vsOthers.bulkRNA != "Inconsistent")
nrow(deg_pbrm1_df)
table(deg_pbrm1_df$PBRM1_Mutated_vsOthers.bulkRNA)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "BAP1_PBRM1_DEGs.logFC.", "BulkRNA.", run_id, ".tsv")
write.table(x = deg_wide_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "BAP1_DEGs.", "BulkRNA.", run_id, ".tsv")
write.table(x = deg_bap1_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "PBRM1_DEGs.", "BulkRNA.", run_id, ".tsv")
write.table(x = deg_pbrm1_df, file = file2write, quote = F, sep = "\t", row.names = F)

