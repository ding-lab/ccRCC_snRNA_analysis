# Yige Wu @ WashU 2019 Oct
## annotate the CPTAC3 ccRCC discovery set cases with mutation status
### mutation: VHL, PBRM1, BAP1, SETD2

# source ------------------------------------------------------------------
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")


# functions ---------------------------------------------------------------
get_yes_no <- function(response_string) {
  simple_reponse <- vector(mode = "character", length = length(response_string)) 
  simple_reponse[grepl(x = response_string, pattern = "Unknown")] <- "Unknown"
  simple_reponse[grepl(x = response_string, pattern = "No")] <- "No"
  simple_reponse[grepl(x = response_string, pattern = "Yes")] <- "Yes"
  return(simple_reponse)
}

# set parameters ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input sample mapping file -----------------------------------------------
clinical_sup_tab <- readxl::read_excel("./Ding_Lab/Projects_Current/CPTAC/CPTAC_Biospecimens_Clinical_Data/CPTAC3_Clinical_Data/CCRCC/CCRCC_Sep2019_flat_file.xlsx")
clinical_sup_tab %>% colnames()

# input the snRNA sample matrix -------------------------------------------
snRNA_sample_mat <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/Sample_Info/03_Post_request_snRNA_sample_selection/RCC_Specimen_Tracking - Case_Matrix.20191007.v1.tsv", data.table = F)

## 65 samples from 30 cases
clinical_tab <- clinical_sup_tab %>%
  filter(case_id %in% snRNA_sample_mat$Case) %>%
  select(-tumor_code) %>%
  rename(Case = case_id) %>%
  rename(Gender = `consent/sex`) %>%
  rename(Age = `consent/age`) %>%
  rename(Race = `consent/ancestry_ethnicity_race_ancestry_identified`) %>%
  rename(Tumor_Stage_Pathological = `baseline/tumor_stage_pathological`) %>%
  rename(BaseLine_Metastasis_Sites = `baseline/specify_distant_metastasis_documented_sites`) %>%
  rename(Baseline_Comments = `baseline/comments`) %>%
  rename(Baseline_IHC = `baseline/ancillary_studies_immunohistochemistry_type_and_results`) %>%
  rename(Adjuvant_Radiation = `follow-up/adjuvant_post_operative_radiation_therapy`) %>%
  rename(Adjuvant_Pharmaceutical_Therapy = `follow-up/adjuvant_post_operative_pharmaceutical_therapy`) %>%
  rename(Adjuvant_immunological_Therapy = `follow-up/adjuvant_post_operative_immunological_therapy`) %>%
  # rename(Cancer_History = `cancer_history/cancer_type`) %>%
  rename(Followup_New_Tumor_Event = `follow-up/new_tumor_event_after_initial_treatment`) %>%
  # mutate(Followup_New_Tumor_Event = get_yes_no(`follow-up/new_tumor_event_after_initial_treatment`)) %>%
  rename(Outcome_of_Initial_Treatment = `follow-up/measure_of_success_of_outcome_at_the_completion_of_initial_first_course_treatment`) %>%
  rename(Outcome_at_Followup_Completion = `follow-up/measure_of_success_of_outcome_at_completion_of_this_follow_up_form`) %>%
  rename(Followup_Cause_of_Death = `follow-up/cause_of_death`) %>%
  rename(Followup_Period = `follow-up/follow_up_period`) %>%
  select(Case, Gender, Age, Race, Tumor_Stage_Pathological, BaseLine_Metastasis_Sites, 
         # Baseline_Comments, Baseline_IHC, 
         Outcome_of_Initial_Treatment, 
         Followup_Period, Outcome_at_Followup_Completion, Followup_Cause_of_Death, Followup_New_Tumor_Event, 
         Adjuvant_Radiation, Adjuvant_Pharmaceutical_Therapy, Adjuvant_immunological_Therapy) %>%
  arrange(Outcome_of_Initial_Treatment, Outcome_at_Followup_Completion)

write.table(x = clinical_tab, file = paste0(dir_out, "snRNA_ccRCC_Clinicl_Table.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)
