# Yige Wu @ WashU 2020 Jul

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
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

# functions ---------------------------------------------------------------
get_yes_no <- function(response_string) {
  simple_reponse <- vector(mode = "character", length = length(response_string)) 
  simple_reponse[grepl(x = response_string, pattern = "Unknown")] <- "Unknown"
  simple_reponse[grepl(x = response_string, pattern = "No")] <- "No"
  simple_reponse[grepl(x = response_string, pattern = "Yes")] <- "Yes"
  return(simple_reponse)
}

# input sample mapping file -----------------------------------------------
clinical_sup_df <- readxl::read_excel("./Resources/Clinical_Data/CCRCC_July2020_clinical_data.xlsx")
colnames_clinical_df <- data.frame(colnames_clinical = colnames(clinical_sup_df))
colnames(clinical_sup_df)
# input the snRNA sample matrix -------------------------------------------
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv", data.table = F)

## rename
clinical_tab <- clinical_sup_df %>%
  filter(case_id %in% idmetadata_df$Case) %>%
  select(-tumor_code) %>%
  rename(Case = case_id) %>%
  rename(Gender = `consent/sex`) %>%
  rename(Age = `consent/age`) %>%
  rename(ethnicity_race_ancestry_identified = `consent/ethnicity_race_ancestry_identified`) %>%
  rename(Race = `consent/race`) %>%
  mutate(Race = ifelse(is.na(Race), ifelse(ethnicity_race_ancestry_identified=="Caucasian", "White", ethnicity_race_ancestry_identified), Race)) %>%
  rename(Ethnicity = `consent/ethnicity`) %>%
  rename(Tumor_Stage_Pathological = `baseline/tumor_stage_pathological`) %>%
  rename(Primary_Tumor_Pathologic_Stage = `baseline/pathologic_staging_primary_tumor`) %>%
  rename(Sarcomatoid_Features = `baseline/sarcomatoid_features`) %>%
  rename(BaseLine_Metastasis_Sites = `baseline/specify_distant_metastasis_documented_sites`) %>%
  rename(Baseline_IHC = `baseline/ancillary_studies_immunohistochemistry_type_and_results`) %>%
  rename(Adjuvant_Radiation = `follow-up/adjuvant_post-operative_radiation_therapy`) %>%
  rename(Adjuvant_Pharmaceutical_Therapy = `follow-up/adjuvant_post-operative_pharmaceutical_therapy`) %>%
  rename(Adjuvant_immunological_Therapy = `follow-up/adjuvant_post-operative_immunological_therapy`) %>%
  # rename(Cancer_History = `cancer_history/cancer_type`) %>%
  rename(Followup_New_Tumor_Event = `follow-up/new_tumor_after_initial_treatment`) %>%
  # mutate(Followup_New_Tumor_Event = get_yes_no(`follow-up/new_tumor_event_after_initial_treatment`)) %>%
  rename(Outcome_of_Initial_Treatment = `follow-up/measure_of_success_of_outcome_at_the_completion_of_initial_first_course_treatment`) %>%
  rename(Outcome_at_Followup_Completion = `follow-up/measure_of_success_of_outcome_at_completion_of_this_follow-up_form`) %>%
  rename(Followup_Cause_of_Death = `follow-up/cause_of_death`) %>%
  rename(Followup_Period = `follow-up/follow_up_period`) %>%
  select(Case, Gender, Age, Race, Ethnicity, ethnicity_race_ancestry_identified,
         Tumor_Stage_Pathological, Primary_Tumor_Pathologic_Stage,
         Sarcomatoid_Features,
         BaseLine_Metastasis_Sites, Baseline_IHC, 
         # Baseline_Comments, 
         Outcome_of_Initial_Treatment, 
         Followup_Period, Outcome_at_Followup_Completion, Followup_Cause_of_Death, Followup_New_Tumor_Event, 
         Adjuvant_Radiation, Adjuvant_Pharmaceutical_Therapy, Adjuvant_immunological_Therapy) %>%
  arrange(Outcome_of_Initial_Treatment, Outcome_at_Followup_Completion)

write.table(x = clinical_tab, file = paste0(dir_out, "snRNA_ccRCC_Clinicl_Table.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)
