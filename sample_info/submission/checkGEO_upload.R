# Yige Wu @ WashU 2021 Nov

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")

# input dependencies ------------------------------------------------------
## input meta data
metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210809.v1/meta_data.20210809.v1.tsv")

# process -----------------------------------------------------------------
aliquots = metadata_df$Aliquot.snRNA[metadata_df$snRNA_available & metadata_df$Case != "C3L-00359"]
aliquots
aliquots_downloaded = list.files(path = "~/Documents/Project/ccRCC_snRNA/Submission/GEO_upload/snRNAseq/")
aliquots_downloaded
aliquots[!(aliquots %in% aliquots_downloaded)]

aliquots = metadata_df$Aliquot.snRNA[metadata_df$snATAC_used & metadata_df$Case != "C3L-00359"]
aliquots
aliquots_downloaded = list.files(path = "~/Documents/Project/ccRCC_snRNA/Submission/GEO_upload/snATACseq/")
aliquots_downloaded
aliquots[!(aliquots %in% aliquots_downloaded)]
