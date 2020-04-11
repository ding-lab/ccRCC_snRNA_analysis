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
celltypemarker_filtered_df <- fread(input = "./Resources/Analysis_Results/findmarkers/findallmarker_roc_subsample_bycelltypeshorter_on_katmai/20200411.v1/findallmarkers_roc_bycelltypeshorter.20200411.v1.tsv", data.table = F)
## input average expression

## input barcode2celltype table

# make adjustment ---------------------------------------------------------
## get unique filtered genes
## filter average expression by filtered genes
## adjust average expression by 

