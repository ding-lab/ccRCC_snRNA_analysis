# Yige Wu @WashU Jul 2020

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

# input bulk & sn profiles ------------------------------------------------
bulk_sn_omicsprofile_df <- fread(input = "./Resources/Analysis_Results/bulk/other/merge_bulk_sn_profiles/20200512.v1/bulk_sn_omics_profile.20200512.v1.tsv", data.table = F)

# input clinical profile --------------------------------------------------
clinical_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/generate_clinical_table/20200717.v1/snRNA_ccRCC_Clinicl_Table.20200717.v1.tsv")


# merge -------------------------------------------------------------------


# write output ------------------------------------------------------------


