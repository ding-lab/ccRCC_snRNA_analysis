# Yige Wu @WashU Dec 2021

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
tpm_df <- fread(input = "./Resources/Analysis_Results/bulk/expression/rna/unite_knockout_cellline_mrna_kallisto/20220510.v1/Knockout_celllines.Kallisto.gene_level.TPM.20220510.v1.tsv", data.table = F)

# process -----------------------------------------------------------------
tpm_caki1_df <- tpm_df[,colnames(tpm_df)[grepl(pattern = "caki1", x = colnames(tpm_df))]]
tpm_caki1_norm_df <- tpm_caki1_df/tpm_caki1_df$caki1_control
tpm_rcc4_df <- tpm_df[,colnames(tpm_df)[grepl(pattern = "rcc4", x = colnames(tpm_df))]]
tpm_rcc4_norm_df <- tpm_rcc4_df/tpm_rcc4_df$rcc4_control
tpm_norm_df <- cbind(tpm_df[, c("gene_name", "ensembl_gene_id")], 
                     tpm_caki1_norm_df[, colnames(tpm_caki1_norm_df)[!grepl(pattern = "control", x = colnames(tpm_caki1_norm_df))]],
                     tpm_rcc4_norm_df[, colnames(tpm_rcc4_norm_df)[!grepl(pattern = "control", x = colnames(tpm_rcc4_norm_df))]])

tpm_norm_filtered_df <- tpm_norm_df %>%
  filter(gene_name %in% c("MXI1", "KLF9", "CP", "PCSK6", "HK2", "PKM", "PFKP", "ENO2", "MYC"))

# write output ------------------------------------------------------------
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
file2write <- paste0(dir_out, "Knockout_celllines.Kallisto.gene_level.TPM.NormalizedByControl.", run_id, ".tsv")
write.table(x = tpm_norm_df, file = file2write, sep = "\t", row.names = F, quote = F)
file2write <- paste0(dir_out, "Knockout_celllines.Kallisto.gene_level.TPM.NormalizedByControl.Filtered.", run_id, ".tsv")
write.table(x = tpm_norm_filtered_df, file = file2write, sep = "\t", row.names = F, quote = F)
