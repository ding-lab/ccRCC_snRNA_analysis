# Yige Wu @WashU Jun 2021

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
## input clinical info
clinical_df <- readxl::read_excel(path =  "~/Box/CPTAC_ccRCC/Data_Freeze_1.1/CPTAC_ccRCC_Combined/Clinical_data/CCRCC_June2021_clinical_data.xlsx")
## input the ccRCC/non-cc RCC info
class_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Case_ID/CPTAC_ccRCC_discovery_caseID_v1.0.tsv")

# extract info ------------------------------------------------------------
## get only ccRCC
case2survival_df <- class_df %>%
  filter(Histologic_Type == "Clear cell renal cell carcinoma") %>%
  select(CASE_ID)
## get the latest survival time
clinical_df$survial_time <- sapply(X = clinical_df$`follow-up/days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact`, function(x) {
  days_vec <- str_split(string = x, pattern = "\\|")[[1]]
  day_last <- days_vec[length(days_vec)]
  return(day_last)
})
clinical_df$with_new_event <- sapply(X = 1:nrow(clinical_df), function(x, days, statuses) {
  status_vec <- str_split(string = statuses[x], pattern = "\\|")[[1]]
  status_last <- status_vec[length(status_vec)]
  
  days_vec <- str_split(string = days[x], pattern = "\\|")[[1]]
  if (!is.na(days_vec) & length(days_vec) >=2) {
    if ((days_vec[length(days_vec)] == days_vec[length(days_vec)-1])) {
      status_last <- status_vec[length(status_vec)-1]
      
    }
  }
  return(status_last)
}, days = clinical_df$`follow-up/days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact`, statuses = clinical_df$`follow-up/tumor_status_at_date_of_last_contact_or_death`)

View(clinical_df[, c("survial_time", "with_new_event",
                     "follow-up/days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact",
                     "follow-up/tumor_status_at_date_of_last_contact_or_death")])

case2survival_df$survival_time <- mapvalues(x = case2survival_df$CASE_ID, from = clinical_df$case_id, to = as.vector(clinical_df$survial_time))
case2survival_df$survival_time <- as.numeric(case2survival_df$survival_time)
case2survival_df$with_new_event <- mapvalues(x = case2survival_df$CASE_ID, from = clinical_df$case_id, to = as.vector(clinical_df$with_new_event))
case2survival_df$with_new_event[case2survival_df$with_new_event == "Unknown"] <- NA

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "CPTAC_Discovery_ccRCC_Survival_Time", run_id, ".tsv")
write.table(file = file2write, x = case2survival_df, quote = F, sep = "\t", row.names = F)
