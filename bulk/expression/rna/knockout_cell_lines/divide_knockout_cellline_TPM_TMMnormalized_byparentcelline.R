# Yige Wu @WashU May 2022

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    print(paste0(pkg_name_tmp, "is being installed!"))
    BiocManager::install(pkgs = pkg_name_tmp, update = F)
    install.packages(pkg_name_tmp, dependencies = T)
  }
  print(paste0(pkg_name_tmp, " is installed!"))
  library(package = pkg_name_tmp, character.only = T)
}


# input dependencies ------------------------------------------------------
exp_df <- fread(input = "./Resources/Analysis_Results/bulk/expression/rna/knockout_cell_lines/normalization/get_edgeR_TMM_normalized_counts/20220511.v1/CPM.TMM_normalized.Knockout_Cell_Lines.20220511.v1.tsv", data.table = F)

# process -----------------------------------------------------------------
exp_caki1_df <- exp_df[,colnames(exp_df)[grepl(pattern = "caki", x = colnames(exp_df))]]
exp_caki1_norm_df <- exp_caki1_df/exp_caki1_df$sample.caki_1_control_e1
exp_rcc4_df <- exp_df[,colnames(exp_df)[grepl(pattern = "rcc4", x = colnames(exp_df))]]
exp_rcc4_norm_df <- exp_rcc4_df/exp_rcc4_df$sample.rcc4_control_e1
# exp_norm_df <- cbind(exp_df[, 1:7], 
#                      exp_caki1_norm_df[, colnames(exp_caki1_norm_df)[!grepl(pattern = "control", x = colnames(exp_caki1_norm_df))]],
#                      exp_rcc4_norm_df[, colnames(exp_rcc4_norm_df)[!grepl(pattern = "control", x = colnames(exp_rcc4_norm_df))]])
exp_norm_df <- cbind(exp_df[, 1:7], 
                     exp_caki1_norm_df,
                     exp_rcc4_norm_df)
exp_norm_filtered_df <- exp_norm_df %>%
  filter(external_gene_name %in% c("MXI1", "KLF9", "CP", "PCSK6", "HK2", "PKM", "PFKP", "ENO2", "MYC"))

# write output ------------------------------------------------------------
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
file2write <- paste0(dir_out, "Knockout_celllines.Kallisto.gene_level.CPM.TMMnormalized.DividedByControl.", run_id, ".tsv")
write.table(x = exp_norm_df, file = file2write, sep = "\t", row.names = F, quote = F)
file2write <- paste0(dir_out, "Knockout_celllines.Kallisto.gene_level.CPM.TMMnormalized.DividedByControl.Filtered.", run_id, ".tsv")
write.table(x = exp_norm_filtered_df, file = file2write, sep = "\t", row.names = F, quote = F)
