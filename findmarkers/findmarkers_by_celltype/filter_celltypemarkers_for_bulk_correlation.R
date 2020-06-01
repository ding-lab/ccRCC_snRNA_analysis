# Yige Wu @WashU Apr 2020

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

# input denpendencies -----------------------------------------------------
## input cell type markers
celltypemarker_df <- fread(input = "./Resources/Analysis_Results/findmarkers/findallmarker_roc_subsample_bycelltypeshorter_on_katmai/20200411.v1/findallmarkers_roc_bycelltypeshorter.20200411.v1.tsv", data.table = F)

# filter ------------------------------------------------------------------
celltypemarker_filtered_df <- celltypemarker_df %>%
  filter(avg_diff > 0) %>%
  filter(power > 0.1)
## the power thershold is an arbitrary cutoff, after looking at IL7R having a predictive power of 0.008 in myofibroblasts
nrow(celltypemarker_filtered_df)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "celltypemarkers_filtered_for_bulk_correlation.", run_id, ".tsv")
write.table(file = file2write, x = celltypemarker_filtered_df, quote = F, sep = "\t", row.names = F)