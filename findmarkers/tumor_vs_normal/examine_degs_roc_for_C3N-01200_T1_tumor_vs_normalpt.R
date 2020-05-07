# Yige Wu @WashU Apr 2020

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

# input dependencies ------------------------------------------------------
## input the hif pathway members
hiftargets_df <- fread(input = "./Resources/Analysis_Results/dependencies/write_hif_targets/20200302.v1/HIF_Target_Genes.20200302.v1.tsv", data.table = F)
## input id meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200413.v1/meta_data.20200413.v1.tsv", data.table = F)
## input DEGs
list.files(path = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/")
degs_df <- fread(input = paste0("./Resources/Analysis_Results/findmarkers/tumor_vs_normal/findmarker_roc_CPT0075130004_tumor_vs_normalpt_on_katmai/20200413.v1/findmarkers_roc.CPT0075130004tumor_vs_normalpt.20200413.v1.tsv"), data.table = F)

# examine genes higher in tumor cells -------------------------------------
degs_higher_in_tumor <- degs_df %>%
  filter(avg_diff > 0) %>%
  filter(row_name %in% hiftargets_df$target_genesymbol)
degs_higher_in_tumor$row_name

degs_lower_in_tumor <- degs_df %>%
  filter(avg_diff < 0) %>%
  filter(row_name %in% hiftargets_df$target_genesymbol)
degs_lower_in_tumor$row_name
