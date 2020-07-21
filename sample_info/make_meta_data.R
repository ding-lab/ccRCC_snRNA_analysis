# Yige Wu @ WashU 2020 March
## make meta data table to map case ID to sample ID to aliquot ID (proteomics + snRNA)
## 2020-04-13 make new aliquot to simplify the aliquot ids

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
## input meta data with proteomics aliquot ids from discovery manuscript
bulk_id_metadata_df <- fread("./Resources/Bulk_Processed_Data/Meta_Data/cptac-metadata.csv")
## input sample info shipped for snRNA at 2019-June
snRNA_id_metadata_df1 <- readxl::read_xlsx("./Resources/Sample_Info/02_Post_request_CPTAC3_shipment/JHU Distribution Manifest_ccRCC material to WUSTL_6.17.2019.xlsx")
## input sample info shipped for snRNA at 2019-Oct
snRNA_id_metadata_df2 <- readxl::read_xlsx("./Resources/Sample_Info/02_Post_request_CPTAC3_shipment/Wash U_ccRCC and GBM bulk tissues_Manifest_10.7.2019.xlsx")
## input info about the original tumor piece for snRNA aliquots
snRNA_id_metadata_df1_original <- readxl::read_xlsx("./Resources/Sample_Info/02_Post_request_CPTAC3_shipment/JHU Distribution Manifest_ccRCC material to WUSTL_6.17.2019_original_segment.xlsx")
## input the new aliquot id for the multi-segment sample
id_multisegment_df <- fread(input = "./Resources/Meta_Data/ccRCC_Specimen_Data_Tracking - Multi-Segment_Info.tsv", data.table = F)
## input samples with snRNA data
snrna_sample_df <- fread(data.table = F, input = "./Resources/snRNA_Processed_Data/scRNA_auto/summary/Seurat_Preprocessing.20200701.v1.tsv")
## input samples with snATAC data
snatac_sample_df <- fread(data.table = F, input = "../ccRCC_snATAC/Resources/snATAC_Data_Generation/snATAC_Sample_Selection_Rationale.tsv")
## input FFPE slides availability
ffpe_sample_df <- readxl::read_excel(path = "./Resources/Sample_Info/02_Post_request_CPTAC3_shipment/Manfiest_ccRCC and GBM slides for IHC @ Wash U_10.28.2019.xlsx")

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

# make new aliquot id -----------------------------------------------------
id_metadata_df <- id_metadata_df %>%
  mutate(Aliquot.snRNA.WU = ifelse(Sample_Type == "Normal", paste0(Case, "-N"),
                                   ifelse(Is_discovery_set, paste0(Case, "-T1"), Aliquot.snRNA)))
id_metadata_df$Aliquot.snRNA.WU <- mapvalues(x = id_metadata_df$Aliquot.snRNA.WU, from = id_multisegment_df$Aliquot, to = as.vector(id_multisegment_df$Naming))
## filter out samples that are not snRNA-seqed
id_metadata_df <- id_metadata_df %>%
  filter(Aliquot.snRNA.WU != Aliquot.snRNA)

# annotate which samples has snRNA data ------------------------
id_metadata_df <- id_metadata_df %>%
  mutate(snRNA_available = (Aliquot.snRNA %in% snrna_sample_df$Aliquot))

# annotate which samples has snATAC data ------------------------
id_metadata_df <- id_metadata_df %>%
  mutate(snATAC_available = (Aliquot.snRNA %in% snatac_sample_df$Aliquot))

# annotate which samples have FFPE slides ---------------------------------
ffpe_sample_df <- as.data.frame(ffpe_sample_df)
ffpe_sample_count_df <- table(ffpe_sample_df$`Corresponding Tissue Block ID`)
ffpe_sample_count_df <- as.data.frame(ffpe_sample_count_df)
id_metadata_df$FFPE_available_slides <- mapvalues(x = id_metadata_df$Sample, from = ffpe_sample_count_df$Var1, to = as.vector(ffpe_sample_count_df$Freq))
id_metadata_df$FFPE_available_slides[id_metadata_df$FFPE_available_slides == id_metadata_df$Sample] <- 0
id_metadata_df$FFPE_available_slides <- as.numeric(id_metadata_df$FFPE_available_slides)

# Write meta data table ---------------------------------------------------
write.table(x = id_metadata_df, file = paste0(dir_out, "meta_data", ".", run_id ,".tsv"), quote = F, sep = "\t", row.names = F)
