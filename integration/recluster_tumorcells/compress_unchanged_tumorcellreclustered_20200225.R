# Yige Wu @WashU Nov 2020

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
## input id meta daa
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")
## input samples whose cell type has been corrected
aliquots_reprocessed_df <- readxl::read_excel(path = "./Resources/snRNA_Processed_Data/Cell_Type_Assignment/Individual_AllClusters/Cells_BySampleByClusterByCellTypeShorter.Over50.20201027.xlsx", sheet = "Sheet1")

# get aliquot ids --------------------------------------------------------------
## path for the 20 corrected tumor samples reclustered 20201109.v1
aliquots_reprocessed <- unique(aliquots_reprocessed_df$aliquot[aliquots_reprocessed_df$Cell_type.shorter.original == "Tumor cells" | aliquots_reprocessed_df$Cell_group.detailed == "Tumor cells"])
aliquots_reprocessed
## path for the  9 uncanged tumor cells reclustered version. 20200225.v1
idmetadata_filtered_df <- idmetadata_df %>%
  filter(Sample_Type == "Tumor")
aliquots_unchanged <- idmetadata_filtered_df$Aliquot.snRNA[!(idmetadata_filtered_df$Aliquot.snRNA %in% aliquots_reprocessed)]
aliquots_unchanged

# input and output --------------------------------------------------------


