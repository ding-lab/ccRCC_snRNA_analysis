# Yige Wu @ WashU 2019 Nov
## make meta data table to map case ID to sample ID to aliquot ID (proteomics + snRNA)

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")


# set run id --------------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)


# set output directory ----------------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# input sample info shipped for snRNA at June-------------------------------------
snRNA_meta_tab1 <- read_xlsx("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Sample_Info/02_Post_request_CPTAC3_shipment/JHU Distribution Manifest_ccRCC material to WUSTL_6.17.2019.xlsx")
snRNA_meta_tab1 <- data.frame(snRNA_meta_tab1)
snRNA_meta_tab1 <- snRNA_meta_tab1 %>%
  select(Subject.ID, Parent.Tumor.Segment.ID, Current.Vial.Label) %>%
  mutate(Aliquot = gsub(x = Current.Vial.Label, pattern = " ", replacement = "")) %>%
  rename(Case = Subject.ID) %>%
  rename(Sample = Parent.Tumor.Segment.ID) %>%
  select(Case, Sample, Aliquot) %>%
  mutate(Sample_Type = "Tumor")

# input sample info shipped for snRNA at Oct-------------------------------------
snRNA_meta_tab2 <- read_xlsx("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Sample_Info/02_Post_request_CPTAC3_shipment/Wash U_ccRCC and GBM bulk tissues_Manifest_10.7.2019.xlsx")
snRNA_meta_tab2 <- data.frame(snRNA_meta_tab2)
snRNA_meta_tab2 <- snRNA_meta_tab2 %>%
  filter(Anatomic.Site == "Kidney") %>%
  mutate(Sample_Type = ifelse(Sample.Type == "NAT", "Normal", "Tumor")) %>%
  mutate(Aliquot = gsub(x = Vial.Label, pattern = " ", replacement = "")) %>%
  select(Subject.ID, Parent.Tissue.ID.s., Aliquot, Sample_Type) %>%
  rename(Case = Subject.ID) %>%
  rename(Sample = Parent.Tissue.ID.s.) %>%
  mutate(Is_discovery_set = T)


# input original segments ------------------------------------------------------------------
snRNA_meta_tab1_original <- read_xlsx("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Sample_Info/02_Post_request_CPTAC3_shipment/JHU Distribution Manifest_ccRCC material to WUSTL_6.17.2019_original_segment.xlsx")
snRNA_meta_tab1 <- snRNA_meta_tab1 %>%
  mutate(Is_discovery_set = (Sample %in% snRNA_meta_tab1_original$`Parent Tumor Segment ID`))
  
# merge the meta tables ---------------------------------------------------
snRNA_meta_tab <- rbind(snRNA_meta_tab1, snRNA_meta_tab2)
snRNA_meta_tab %>% head()

# input meta data from discovery manuscript -------------------------------
bulk_meta_tab <- fread("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/cptac-metadata.csv")
bulk_meta_tab %>% head()
bulk_meta_tab <- bulk_meta_tab %>%
  mutate(Is_discovery_set = (Set.A == "yes")) %>%
  select(Case.ID, Specimen.Label, Type, Is_discovery_set) %>%
  rename(Case = Case.ID) %>%
  rename(Aliquot = Specimen.Label) %>%
  rename(Sample_Type =  Type)

# merge snRNA meta data with proteomics meta data -------------------------
meta_tab <- merge(snRNA_meta_tab, bulk_meta_tab, by = c("Case", "Sample_Type", "Is_discovery_set"), all.x = T, suffixes = c(".snRNA", ".bulk"))

# Write meta data table ---------------------------------------------------
write.table(x = meta_tab, file = paste0(dir_out, "meta_data", ".", run_id ,".tsv"), quote = F, sep = "\t", row.names = F)
