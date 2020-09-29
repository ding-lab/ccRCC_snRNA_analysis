# Yige Wu @WashU Sep 2020
## reference of the cellphonedb output: https://www.cellphonedb.org/documentation

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
## input cellphonedb output
cellphone_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/summarize_interactions/summarize_across_samples_by_pair/20200923.v1/cellphonedb.summary_across_samples_by_pair.20200923.v1.tsv")

# filter ------------------------------------------------------------------
colnames(cellphone_df)
cellphone_filtered_df <- cellphone_df %>%
  filter(paired_cellgroups.general == "Stroma&Tumor") %>%
  filter(number_sig_cases >= 20) %>%
  arrange(desc(number_sig_cases), avg_rank_sig_mean.paired_cellgroups.general)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "cellphonedb.summary_across_samples_by_pair.", "stroma&tumor.filtered.", run_id, ".tsv")
write.table(x = cellphone_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
