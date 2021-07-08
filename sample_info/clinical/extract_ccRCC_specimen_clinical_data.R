# Yige Wu @ WashU 2020 Jul
## annotate the CPTAC3 ccRCC discovery set cases with mutation status

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

# functions ---------------------------------------------------------------
get_yes_no <- function(response_string) {
  simple_reponse <- vector(mode = "character", length = length(response_string)) 
  simple_reponse[grepl(x = response_string, pattern = "Unknown")] <- "Unknown"
  simple_reponse[grepl(x = response_string, pattern = "No")] <- "No"
  simple_reponse[grepl(x = response_string, pattern = "Yes")] <- "Yes"
  return(simple_reponse)
}

# input sample mapping file -----------------------------------------------
clinical_sup_df <- readxl::read_excel("./Resources/Clinical_Data/CCRCC_June2021_clinical_data.xlsx", sheet = "Specimen")
colnames_clinical_df <- data.frame(colnames_clinical = colnames(clinical_sup_df))
colnames(clinical_sup_df)
# input the snRNA sample matrix -------------------------------------------
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210423.v1/meta_data.20210423.v1.tsv", data.table = F)

# merge -------------------------------------------------------------------
clinical_df <- merge(idmetadata_df, clinical_sup_df, by.x = c("Case", "Sample"), by.y = c("case_id", "specimen_id"), all = T)

# mutate and rename -------------------------------------------------------
clinical_df$`cptac_path/histologic_grade`
clinical_selected_df <- clinical_df %>%
  mutate(Histologic_Grade_CPTAC = str_split_fixed(string = `cptac_path/histologic_grade`, pattern = ": ", n = 2)[,1]) %>%
  mutate(Histologic_Grade = Histologic_Grade_CPTAC) %>%
  mutate(Histologic_Grade_Description = str_split_fixed(string = `cptac_path/histologic_grade`, pattern = ": ", n = 2)[,2]) %>%
  rename(Histologic_Type_Of_Normal_Tissue = `cptac_path/histologic_type_of_normal_tissue`) %>%
  rename(Histologic_Type = `cptac_path/histologic_type`) %>%
  mutate(Histologic_Type = ifelse(!is.na(Histologic_Type), Histologic_Type, "Normal Adjacent Tissue")) %>%
  rename(Tissue_Type = `tissue_segment/tissue_type`) %>%
  select(Case, Sample, Sample_Type, Tissue_Type, Aliquot.snRNA, Aliquot.bulk, Aliquot.snRNA.WU, discovery_study, confirmatory_study,
         snRNA_available, snATAC_available,
         Histologic_Type, Histologic_Type_Of_Normal_Tissue, Histologic_Grade_CPTAC,
         Histologic_Grade, Histologic_Grade_Description)

# correct -----------------------------------------------------------------
## correct histological grade with Siqi's evaluation
clinical_selected_df$Histologic_Grade[clinical_selected_df$Sample == "C3N-01200-02"] <- "G2"
clinical_selected_df$Histologic_Grade[clinical_selected_df$Sample == "C3N-01200-03"] <- "G4"
clinical_selected_df$Histologic_Grade[clinical_selected_df$Sample == "C3N-01200-01"] <- "G4"
## correct normal samples
clinical_selected_df$Histologic_Grade[clinical_selected_df$Tissue_Type %in% c("normal", "blood")] <- "NAT"
# clinical_selected_df$Histologic_Grade
# clinical_selected_df$Histologic_Type

# write output ------------------------------------------------------------
write.table(x = clinical_selected_df, file = paste0(dir_out, "ccRCC_Specimen_Clinicl_Data.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)
