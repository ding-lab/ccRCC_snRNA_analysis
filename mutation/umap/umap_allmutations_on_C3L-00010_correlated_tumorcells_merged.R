# Yige Wu @WashU Aug 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies-----------------------------------------
## input umap data
umap_df <- fread(data.table = F, input = "./Resources/Analysis_Results/fetch_data/fetchdata_tumorcells_C3L-00010_C3L-00416_C3L-00583_C3N-00242_merged_on_katmai/20200817.v1/C3L-00010_C3L-00416_C3L-00583_C3N-00242.Metadata.20200817.v1.tsv")
## input metadata
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")
