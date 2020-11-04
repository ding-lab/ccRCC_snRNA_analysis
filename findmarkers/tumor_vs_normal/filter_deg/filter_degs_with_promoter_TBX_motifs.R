# Yige Wu @WashU Oct 2020

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
## input the DEG-TF matrix
deg2tf_wide_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/annotate_deg/annotate_degs_associated_with_tf_with_bulk_expression/20201027.v1/DEGs_with_TFs_inDARs.sn_and_Bulk_DE_Annotated.tsv")

# filter to those with at least one TBX motif in promoter -----------------
deg2tf_filtered_df <- deg2tf_wide_df %>%
  filter(MGA == 1 | TBX1 == 1 | TBX15 == 1 | TBX3 == 1 | TBX4 == 1 | TBX5 == 1 | TBX6 == 1)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "DEGs_with_TBXs_inDARs", ".Promoter", ".sn_and_Bulk_DE_Annotated", ".tsv")
write.table(x = deg2tf_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)

