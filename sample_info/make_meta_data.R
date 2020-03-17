# Yige Wu @ WashU 2020 March
## make meta data table to map case ID to sample ID to aliquot ID (proteomics + snRNA)
## 2020-03-16 make new aliquot to simplify the aliquot ids

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
library(dplyr)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input meta data with proteomics aliquot ids from discovery manuscript
bulk_id_metadata_df <- fread("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/cptac-metadata.csv")
## input sample info shipped for snRNA at 2019-June
snRNA_id_metadata_df1 <- read_xlsx("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Sample_Info/02_Post_request_CPTAC3_shipment/JHU Distribution Manifest_ccRCC material to WUSTL_6.17.2019.xlsx")
## input sample info shipped for snRNA at 2019-Oct
snRNA_id_metadata_df2 <- read_xlsx("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Sample_Info/02_Post_request_CPTAC3_shipment/Wash U_ccRCC and GBM bulk tissues_Manifest_10.7.2019.xlsx")
## input info about the original tumor piece for snRNA aliquots
snRNA_id_metadata_df1_original <- read_xlsx("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Sample_Info/02_Post_request_CPTAC3_shipment/JHU Distribution Manifest_ccRCC material to WUSTL_6.17.2019_original_segment.xlsx")

# transform snRNA aliquot id table-------------------------------------
## transform sample info shipped for snRNA at 2019-June
snRNA_id_metadata_df1 <- data.frame(snRNA_id_metadata_df1)
snRNA_id_metadata_df1 <- snRNA_id_metadata_df1 %>%
  dplyr::select(Subject.ID, Parent.Tumor.Segment.ID, Current.Vial.Label) %>%
  dplyr::mutate(Aliquot = gsub(x = Current.Vial.Label, pattern = " ", replacement = "")) %>%
  dplyr::rename(Case = Subject.ID) %>%
  dplyr::rename(Sample = Parent.Tumor.Segment.ID) %>%
  dplyr::select(Case, Sample, Aliquot) %>%
  dplyr::mutate(Sample_Type = "Tumor") %>%
  mutate(Is_discovery_set = (Sample %in% snRNA_id_metadata_df1_original$`Parent Tumor Segment ID`))

## transform sample info shipped for snRNA at 2019-Oct
snRNA_id_metadata_df2 <- data.frame(snRNA_id_metadata_df2)
snRNA_id_metadata_df2 <- snRNA_id_metadata_df2 %>%
  dplyr::filter(Anatomic.Site == "Kidney") %>%
  dplyr::mutate(Sample_Type = ifelse(Sample.Type == "NAT", "Normal", "Tumor")) %>%
  dplyr::mutate(Aliquot = gsub(x = Vial.Label, pattern = " ", replacement = "")) %>%
  dplyr::select(Subject.ID, Parent.Tissue.ID.s., Aliquot, Sample_Type) %>%
  dplyr::rename(Case = Subject.ID) %>%
  dplyr::rename(Sample = Parent.Tissue.ID.s.) %>%
  dplyr::mutate(Is_discovery_set = T)
## merge the meta tables
snRNA_id_metadata_df <- rbind(snRNA_id_metadata_df1, snRNA_id_metadata_df2)
snRNA_id_metadata_df %>% head()

# merge snRNA meta data with proteomics meta data -------------------------
## filter the proteomics aliquot ids
bulk_id_metadata_df <- bulk_id_metadata_df %>%
  mutate(Is_discovery_set = (Set.A == "yes")) %>%
  select(Case.ID, Specimen.Label, Type, Is_discovery_set) %>%
  rename(Case = Case.ID) %>%
  rename(Aliquot = Specimen.Label) %>%
  rename(Sample_Type =  Type)
## merge
id_metadata_df <- merge(snRNA_id_metadata_df, bulk_id_metadata_df, by = c("Case", "Sample_Type", "Is_discovery_set"), all.x = T, suffixes = c(".snRNA", ".bulk"))

# Write meta data table ---------------------------------------------------
write.table(x = id_metadata_df, file = paste0(dir_out, "meta_data", ".", run_id ,".tsv"), quote = F, sep = "\t", row.names = F)
