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
exp_df <- fread("./Resources/Bulk_Processed_Data/Protein/6_CPTAC3_CCRCC_Whole_abundance_gene_protNorm=2_CB.tsv", data.table = F)
## input bulk meta data
metadata_bulk_df <- fread("./Resources/Bulk_Processed_Data/Meta_Data/cptac-metadata.csv")
## input clinical info
case_clinical_df <- readxl::read_excel(path = "./Resources/Bulk_Processed_Data/CPTAC3-ccRCC-SupplementaryTables_Final/Table S1.xlsx", sheet = "ccrcc_clinical_characteristics")

# get the aliquot IDs for bulk corresponding to the snRNA aliquots --------
metadata_filtered_df <- metadata_bulk_df %>%
  filter(Set.A == "yes") %>%
  filter(Specimen.Label != "CPT0012090003")
metadata_filtered_df$Histologic_Type <- mapvalues(x = metadata_filtered_df$Case.ID, from = case_clinical_df$Case_ID, to = case_clinical_df$Histologic_Type)
metadata_filtered_df <- metadata_filtered_df %>%
  filter(Histologic_Type == "Clear cell renal cell carcinoma")
ids_exp_normal <- metadata_filtered_df$Specimen.Label[metadata_filtered_df$Type == "Normal"]
ids_exp_normal

ids_case <- mapvalues(x = ids_exp_normal, from = metadata_filtered_df$Specimen.Label, to = metadata_filtered_df$Case.ID)
ids_case

ids_exp_tumor <- mapvalues(x = ids_case, from = metadata_filtered_df$Case.ID[metadata_filtered_df$Type == "Tumor"], to = as.vector(metadata_filtered_df$Specimen.Label[metadata_filtered_df$Type == "Tumor"]))
ids_exp_tumor

# average for tumors ------------------------------------------------------
exp_tumor_df <- 


