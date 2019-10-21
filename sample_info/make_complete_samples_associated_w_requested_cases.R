# Yige Wu @ WashU 2019 Aug
## make the files per Li's request for sample request:
# Yige,
# 
# Please send me the list. Two files: 1. complete (all samples related to a case) with to be requested highlighted; 2. Just what we are requesting.
# 
# Thanks
# 
# Li

# source ------------------------------------------------------------------
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")


# input CPTAC specimen info -----------------------------------------------
cptac_ccRCC_specimen_tab <- readxl::read_excel(path = "./Ding_Lab/Projects_Current/CPTAC/CPTAC_Biospecimens_Clinical_Data/CPTAC3_Clinical_Data/CCRCC/CCRCC_Sep2019_flat_file.xlsx", sheet = "Specimen")
cptac_ccRCC_specimen_tab %>% colnames()

# input requested sample matrix -------------------------------------------
sample_request_mat <- readxl::read_excel(path = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/Sample_Info/00_Pre_request_sample_selection/WashU_ccRCC_Request_Sample_Matrix.092319.v1.xlsx")
discovery_tumor_request_case <- sample_request_mat$Case[grepl(pattern = "To-be-requested", x = sample_request_mat$Discovery_set_Tumor_segment, ignore.case = T)]
discovery_tumor_request_case

normal_request_case <- sample_request_mat$Case[grepl(pattern = "To-be-requested", x = sample_request_mat$NAT, ignore.case = T)]
normal_request_case

discovery_normal_request_case <- normal_request_case[normal_request_case %in% cptac_ccRCC_specimen_tab$case_id[cptac_ccRCC_specimen_tab$discovery_study == "Yes" & cptac_ccRCC_specimen_tab$`tissue_segment/tissue_type` == "normal"]]
confirmatory_normal_request_case <- normal_request_case[!(normal_request_case %in% discovery_normal_request_case)]
confirmatory_normal_request_case


# filter CPTAC inventory by requested cases ----------------------------------------------------------
requested_cptac_ccRCC_specimen_tab <- cptac_ccRCC_specimen_tab %>%
  filter(case_id %in% sample_request_mat$Case) %>%
  mutate(washu_requested = (`tissue_segment/tissue_type` == "tumor" & discovery_study == "Yes" & case_id %in% discovery_tumor_request_case) | (`tissue_segment/tissue_type` == "normal" & discovery_study == "Yes" & case_id %in% discovery_normal_request_case) | (`tissue_segment/tissue_type` == "normal" & discovery_study == "No" & case_id %in% confirmatory_normal_request_case)) %>%
  select(tumor_code, case_id, specimen_id, washu_requested, qc_cassette_id, slide_id, discovery_study, `tissue_segment/tissue_type`, `tissue_segment/segment`, `tissue_segment/weight_mg`)
write.table(x = requested_cptac_ccRCC_specimen_tab, file = paste0(makeOutDir(), "WashU_ccRCC_Request_Case_List.092319.v1.tsv"), quote = F, sep = "\t", row.names = F)
