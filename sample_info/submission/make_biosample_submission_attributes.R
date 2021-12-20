# Yige Wu @ WashU 2021 Nov

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
## input meta data
metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210809.v1/meta_data.20210809.v1.tsv")
## input clinical data
clinical_case_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/clinical/extract_case_clinical_data_snrna_samples/20211011.v1/snRNA_ccRCC_Clinicl_Table.20211011.v1.tsv")

# process -----------------------------------------------------------------
metadata_filtered_df <- metadata_df %>%
  filter(Case %in% clinical_case_df$Case) %>%
  filter(snRNA_available) %>%
  select(Aliquot.snRNA.WU, Case, snRNA_available, snATAC_used, Sample_Type)
attributes_df <- rbind(metadata_filtered_df %>%
                         mutate(sample_name = paste0("snRNA.", Aliquot.snRNA.WU)) %>%
                         select(sample_name, Case, Sample_Type),
                       metadata_filtered_df %>%
                         filter(snATAC_used) %>%
                         mutate(sample_name = paste0("snATAC.", Aliquot.snRNA.WU)) %>%
                         select(sample_name, Case, Sample_Type))
attributes_df <- merge(x = attributes_df, y = clinical_case_df, by = c("Case"), all.x = T)
attributes_df <- attributes_df %>%
  mutate(sample_title = "") %>%
  mutate(bioproject_accession = "") %>%
  mutate(organism = "Homo sapiens") %>%
  mutate(isolate = Case) %>%
  rename(age = Age) %>%
  mutate(sex = tolower(Gender)) %>%
  mutate(biomaterial_provider = "CPTAC Biospecimen Core Resource, Van Andel Research Institute, 333 Bostwick Ave, Grand Rapids MI 49503") %>%
  mutate(tissue = ifelse(Sample_Type == "Tumor", "primary_tumor", "normal_adjacent_tissue")) %>%
  mutate(cell_line = "") %>%
  mutate(cell_subtype = "") %>%
  mutate(cell_type = "") %>%
  mutate(culture_collection = "") %>%
  mutate(dev_stage = "") %>%
  mutate(disease = "clear cell renal cell carcinoma") %>%
  mutate(disease_stage = "Primary") %>%
  mutate(ethnicity = "") %>%
  mutate(health_state = "") %>%
  mutate(karyotype = "") %>%
  mutate(phenotype = "") %>%
  mutate(population = "") %>%
  rename(race = Race) %>%
  mutate(sample_type = "") %>%
  mutate(treatment = "") %>%
  mutate(description = "") %>%
  select(sample_name, sample_title, bioproject_accession, organism, isolate, age, biomaterial_provider, sex, tissue,
         cell_line, cell_subtype, cell_type, culture_collection, dev_stage, disease, disease_stage,
         ethnicity, health_state, karyotype, phenotype, population, race, sample_type, treatment, description)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "HumanTissue.SubmissionAttributes.", run_id, ".tsv")
write.table(x = attributes_df, file = file2write, quote = F, sep = "\t", row.names = F)


