# Yige Wu @WashU Mar 2022
## reference: https://www.emilyzabor.com/tutorials/survival_analysis_in_r_tutorial.html

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}

# input dependencies ------------------------------------------------------
## input clinical info
clinical_df <- readxl::read_excel(path =  "../../../../CPTAC_ccRCC/Data_Freeze_1.1/CPTAC_ccRCC_Combined/Clinical_data/CCRCC_Oct2021_clinical_data.xlsx")
## input the ccRCC/non-cc RCC info
class_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Case_ID/CPTAC_ccRCC_discovery_caseID_v1.0.tsv")
## inpu the specimen clinical info
clinical_specimen_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/clinical/extract_ccRCC_specimen_clinical_data/20210706.v1/ccRCC_Specimen_Clinicl_Data.20210706.v1.tsv")

# extract info ------------------------------------------------------------
colnames(clinical_df)
clinical_df$`follow-up/vital_status_at_date_of_last_contact`
clinical_df$`follow-up/days_from_date_of_initial_pathologic_diagnosis_to_date_of_death`
clinical_df$`follow-up/new_tumor_after_initial_treatment`
clinical_df$`follow-up/days_from_date_of_initial_pathologic_diagnosis_to_date_of_new_tumor_after_initial_treatment`
clinical_df$`follow-up/measure_of_success_of_outcome_at_date_of_last_contact_or_death`
clinical_df$`follow-up/measure_of_success_of_outcome_at_the_completion_of_initial_first_course_treatment`

## get the latest survival time
clinical_df$days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact <- sapply(X = clinical_df$`follow-up/days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact`, function(x) {
  days_vec <- str_split(string = x, pattern = "\\|")[[1]]
  day_last <- days_vec[length(days_vec)]
  return(day_last)
})
clinical_df$vital_status_at_date_of_last_contact <- sapply(X = clinical_df$`follow-up/vital_status_at_date_of_last_contact`, function(x) {
  tmp_vec <- str_split(string = x, pattern = "\\|")[[1]]
  tmp <- tmp_vec[length(tmp_vec)]
  return(tmp)
})
clinical_df$days_from_date_of_initial_pathologic_diagnosis_to_date_of_death <- sapply(X = clinical_df$`follow-up/days_from_date_of_initial_pathologic_diagnosis_to_date_of_death`, function(x) {
  tmp_vec <- str_split(string = x, pattern = "\\|")[[1]]
  tmp <- tmp_vec[length(tmp_vec)]
  return(tmp)
})
clinical_df$new_tumor_after_initial_treatment <- sapply(X = clinical_df$`follow-up/new_tumor_after_initial_treatment`, function(x) {
  tmp_vec <- str_split(string = x, pattern = "\\|")[[1]]
  tmp <- tmp_vec[length(tmp_vec)]
  return(tmp)
})
clinical_df$days_from_date_of_initial_pathologic_diagnosis_to_date_of_new_tumor_after_initial_treatment <- sapply(X = clinical_df$`follow-up/days_from_date_of_initial_pathologic_diagnosis_to_date_of_new_tumor_after_initial_treatment`, function(x) {
  tmp_vec <- str_split(string = x, pattern = "\\|")[[1]]
  tmp <- tmp_vec[length(tmp_vec)]
  return(tmp)
})
clinical_df$tumor_status_at_date_of_last_contact_or_death <- sapply(X = clinical_df$`follow-up/tumor_status_at_date_of_last_contact_or_death`, function(x) {
  tmp_vec <- str_split(string = x, pattern = "\\|")[[1]]
  tmp <- tmp_vec[length(tmp_vec)]
  return(tmp)
})
clinical_df$days_from_date_of_initial_pathologic_diagnosis_to_withtumor <- sapply(X = 1:nrow(clinical_df), function(i, fu_time, tumor_status) {
  tumorstatus_vec <- str_split(string = tumor_status[i], pattern = "\\|")[[1]]
  fudays_vec <- str_split(string = fu_time[i], pattern = "\\|")[[1]]
  if (!is.na(tumorstatus_vec)  & any(tumorstatus_vec == "With Tumor")) {
    tmp <- fudays_vec[which(tumorstatus_vec == "With Tumor")[1]]
  } else {
    tmp <- NA
  }
  return(tmp)
}, 
fu_time = clinical_df$`follow-up/days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact`, 
tumor_status = clinical_df$`follow-up/tumor_status_at_date_of_last_contact_or_death`)

as.vector(clinical_df$days_from_date_of_initial_pathologic_diagnosis_to_withtumor)

# clinical_df %>%
#   select(days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact, 
#          vital_status_at_date_of_last_contact, 
#          days_from_date_of_initial_pathologic_diagnosis_to_date_of_death,
#          new_tumor_after_initial_treatment,
#          days_from_date_of_initial_pathologic_diagnosis_to_date_of_new_tumor_after_initial_treatment,
#          tumor_status_at_date_of_last_contact_or_death,
#          `follow-up/days_from_date_of_initial_pathologic_diagnosis_to_date_of_additional_surgery_for_new_tumor_metastasis`) %>%
#   View()
## there are some cases where there was new tumor after initial treatment, but at the last contact or death, they are tumor free
clinical_df %>%
  filter(tumor_status_at_date_of_last_contact_or_death == "With Tumor" & new_tumor_after_initial_treatment == "No") %>%
  View()

clinical_df %>%
  filter(tumor_status_at_date_of_last_contact_or_death == "With Tumor") %>%
  View()

# make data frame ---------------------------------------------------------
## add basic information for future survival analysis
### get only ccRCC
case2survival_df <- class_df %>%
  filter(Histologic_Type == "Clear cell renal cell carcinoma") %>%
  select(CASE_ID)
case2survival_df$age <- mapvalues(x = case2survival_df$CASE_ID, from = clinical_df$case_id, to = as.vector(clinical_df$`consent/age`))
case2survival_df$age[case2survival_df$age == ">=90"] <- "90"
case2survival_df$age <- as.numeric(case2survival_df$age)
case2survival_df$sex <- mapvalues(x = case2survival_df$CASE_ID, from = clinical_df$case_id, to = as.vector(clinical_df$`consent/sex`))
case2survival_df$sex.ismale <- as.numeric(case2survival_df$sex == "Male")
case2survival_df$tumor_stage_pathological <- mapvalues(x = case2survival_df$CASE_ID, from = clinical_df$case_id, to = as.vector(clinical_df$`baseline/tumor_stage_pathological`))
case2survival_df <- case2survival_df %>%
  mutate(tumor_stage_pathological = gsub(x = tumor_stage_pathological, pattern = "Stage ", replacement = "")) %>%
  mutate(stage.numeric = ifelse(tumor_stage_pathological == "I", 1,
                                ifelse(tumor_stage_pathological == "II", 2,
                                       ifelse(tumor_stage_pathological == "III", 3, 4))))
### add tumor grade
clinical_specimen_filtered_df <- clinical_specimen_df %>%
  filter(Tissue_Type == "tumor") %>%
  mutate(Histologic_Grade.numeric = gsub(x = Histologic_Grade, pattern = "G", replacement = "")) %>%
  arrange(Case, desc(Histologic_Grade.numeric)) %>%
  select(Case, Histologic_Grade.numeric)
clinical_specimen_filtered_df <- clinical_specimen_filtered_df[!duplicated(clinical_specimen_filtered_df$Case),]
case2survival_df$grade.numeric <- mapvalues(x = case2survival_df$CASE_ID, from = clinical_specimen_filtered_df$Case, to = as.vector(clinical_specimen_filtered_df$Histologic_Grade.numeric))
case2survival_df$grade.numeric <- as.numeric(case2survival_df$grade.numeric)
## add survival time
case2survival_df <- merge(x = case2survival_df, y = clinical_df %>%
                            select(case_id, days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact, 
                                   vital_status_at_date_of_last_contact,
                                   days_from_date_of_initial_pathologic_diagnosis_to_date_of_death,
                                   days_from_date_of_initial_pathologic_diagnosis_to_withtumor,
                                   new_tumor_after_initial_treatment,
                                   days_from_date_of_initial_pathologic_diagnosis_to_date_of_new_tumor_after_initial_treatment,
                                   tumor_status_at_date_of_last_contact_or_death),
                          by.x = c("CASE_ID"), by.y = c("case_id"), all.x = T)
case2survival_df$days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact <- as.numeric(case2survival_df$days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact)
case2survival_df$days_from_date_of_initial_pathologic_diagnosis_to_date_of_death <- as.numeric(case2survival_df$days_from_date_of_initial_pathologic_diagnosis_to_date_of_death)
case2survival_df$days_from_date_of_initial_pathologic_diagnosis_to_withtumor <- as.numeric(case2survival_df$days_from_date_of_initial_pathologic_diagnosis_to_withtumor)

case2survival_df$tumor_status_at_date_of_last_contact_or_death

# make components for overall survival ------------------------------------
case2survival_df <- case2survival_df %>%
  mutate(OS_time = days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact) %>%
  mutate(OS_status = ifelse(vital_status_at_date_of_last_contact == "Deceased", "dead", "censored")) %>%
  mutate(PFS_time = ifelse(!is.na(days_from_date_of_initial_pathologic_diagnosis_to_withtumor), days_from_date_of_initial_pathologic_diagnosis_to_withtumor,
                           days_from_date_of_initial_pathologic_diagnosis_to_date_of_last_contact)) %>%
  mutate(PFS_status = ifelse(!is.na(days_from_date_of_initial_pathologic_diagnosis_to_withtumor) | vital_status_at_date_of_last_contact == "Deseased", "progression", "censored"))


# write output ------------------------------------------------------------
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
file2write <- paste0(dir_out, "CPTAC_Discovery_ccRCC_Survival_Time", run_id, ".tsv")
write.table(file = file2write, x = case2survival_df, quote = F, sep = "\t", row.names = F)
