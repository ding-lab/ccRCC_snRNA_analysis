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
deg_snRNA_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/pbrm1_bap1_vs_non_mutants/summarize_degs/summarize_PBRM1_BAP1_DEGs/20210405.v1/BAP1_DEGs.20210405.v1.tsv")
deg_bulkRNA_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/edgeR/examine_degs/summarize_PBRM1_BAP1_bulkRNA_DEGs/20210405.v1/BAP1_PBRM1_DEGs.logFC.BulkRNA.20210405.v1.tsv")

# overlap -----------------------------------------------------------------
deg_merged_df <- merge(x = deg_snRNA_df, y = deg_bulkRNA_df, by = c("genesymbol_deg"), all.x = T)
deg_merged_filtered_df <- deg_merged_df %>%
  filter(!is.na(BAP1_Mutated_vsOthers.bulkRNA) & BAP1_vs_PBRM1_Mutants_snRNA == BAP1_Mutated_vsOthers.bulkRNA)
deg_merged_filtered_df <- deg_merged_filtered_df %>%
  arrange(desc(BAP1_vs_PBRM1_Mutants_Tumorcells_Up), desc(BAP1_vs_NonMutants_Tumorcells_Up), desc(BAP1_Tumorcells_vs_PTcells_Up))
table(deg_merged_filtered_df$BAP1_vs_PBRM1_Mutants_snRNA)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "BAP1_DEGs.", "snRNA_and_bulkRNA.", run_id, ".tsv")
write.table(x = deg_merged_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "BAP1_DEGs.Overlap.", "snRNA_and_bulkRNA.", run_id, ".tsv")
write.table(x = deg_merged_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)