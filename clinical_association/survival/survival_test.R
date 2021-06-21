# Yige Wu @WashU Jun 2021

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
## input survival time
survival_df <- readxl::read_excel(path =  "~/Box/CPTAC_ccRCC/Data_Freeze_1.1/CPTAC_ccRCC_Combined/Clinical_data/CCRCC_Apr2021_clinical_data.xlsx")
survival_df <- readxl::read_excel(path =  "~/Box/CPTAC_ccRCC/Data_Freeze_1.1/CPTAC_ccRCC_Discovery/Case_ID/CPTAC_ccRCC_discovery_caseID_v1.0.tsv")

survival_df$`follow-up/days_from_date_of_initial_pathologic_diagnosis_to_date_of_death`
survival_df$`follow-up/days_from_date_of_collection_to_date_of_last_contact`

View(survival_df[, c("follow-up/days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact",
                     "follow-up/days_from_date_of_collection_to_date_of_last_contact", 
                     "follow-up/days_from_date_of_initial_pathologic_diagnosis_to_date_of_death")])

View(survival_df[, c("follow-up/days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact",
                     "follow-up/tumor_status_at_date_of_last_contact_or_death", 
                     "follow-up/measure_of_success_of_outcome_at_date_of_last_contact_or_death" )])
