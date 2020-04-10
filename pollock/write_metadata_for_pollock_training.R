# Yige Wu @WashU March 2020
## for each individual sample, write cell meta data

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the cell to cell type table
barcode2celltype_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/map_celltype_to_barcode/20200224.v1/30_aliquot_integration.barcode2celltype.20200224.v1.tsv", data.table = F)

# add more info -----------------------------------------------------------
## add more info
pollock_metadata_df <- barcode2celltype_df %>%
  rename(cell_id = individual_barcode) %>%
  rename(sample_id = orig.ident) %>%
  mutate(cancer_type = "KIRC") %>%
  mutate(tissue_type = "kidney") %>%
  mutate(organ_type = "kidney") %>%
  mutate(species = "Homo sapiens") %>%
  mutate(method = "sn") %>%
  mutate(facs = "yes") %>%
  mutate(cell_type = Most_Enriched_Cell_Type4)
## make cell type annotation
pollock_metadata_df$cell_type[pollock_metadata_df$cell_type == ""] <- pollock_metadata_df$Most_Enriched_Cell_Type3[pollock_metadata_df$cell_type == ""]
pollock_metadata_df$cell_type[pollock_metadata_df$cell_type == ""] <- pollock_metadata_df$Most_Enriched_Cell_Type2[pollock_metadata_df$cell_type == ""]
pollock_metadata_df$cell_type[pollock_metadata_df$cell_type == ""] <- pollock_metadata_df$Most_Enriched_Cell_Type1[pollock_metadata_df$cell_type == ""]
pollock_metadata_df$cell_type[pollock_metadata_df$cell_type == ""] <- pollock_metadata_df$Most_Enriched_Cell_Group[pollock_metadata_df$cell_type == ""]
table(pollock_metadata_df$cell_type[pollock_metadata_df$Is_Normal_Nephron_Epithelium])

## correct cell type annotation according to:
## https://docs.google.com/spreadsheets/d/10qV-4VoaVPRfd5HdAqsGo5omv_ImlOc011pvkP7GqH0/edit#gid=712436422
pollock_metadata_df$cell_type[!pollock_metadata_df$Is_Normal_Nephron_Epithelium & pollock_metadata_df$Most_Enriched_Cell_Group == "Nephron_Epithelium"] <- "Malignant proximal tubule"
pollock_metadata_df$cell_type[pollock_metadata_df$cell_type == "cDC1"] <- "DC"
pollock_metadata_df$cell_type[pollock_metadata_df$cell_type == "Myofibroblasts"] <- "Fibroblasts"
pollock_metadata_df$cell_type[pollock_metadata_df$cell_type == "NK-cells"] <- "NK cells"
pollock_metadata_df$cell_type[pollock_metadata_df$cell_type == "Thick Ascending Limb"] <- "Ascending Loop of Henle"
table(pollock_metadata_df$cell_type)

## select columns
## according to https://docs.google.com/spreadsheets/d/10qV-4VoaVPRfd5HdAqsGo5omv_ImlOc011pvkP7GqH0/edit#gid=1547742495
pollock_metadata_df <- pollock_metadata_df %>%
  select(cell_id, sample_id, cancer_type, tissue_type, organ_type, cell_type, species, method, facs)

# write table -------------------------------------------------------------
write.table(x = pollock_metadata_df, file = paste0(dir_out, "ccRCC_snRNA_Pollock_Cell_Metadata.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)

