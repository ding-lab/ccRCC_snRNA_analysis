# Yige Wu @ WashU 2019 Sep
## make meta data table to map case ID to sample ID to aliquot ID (proteomics + snRNA)


# source ------------------------------------------------------------------
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")


# input sample info shipped for snRNA -------------------------------------
snRNA_meta_tab <- read_xlsx("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/Sample_Info/02_Post_request_CPTAC3_shipment/JHU Distribution Manifest_ccRCC material to WUSTL_6.17.2019_original_segment.xlsx")
snRNA_meta_tab <- data.frame(snRNA_meta_tab)
snRNA_meta_tab %>% head()
# input meta data from discovery manuscript -------------------------------
bulk_meta_tab <- fread("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/cptac-metadata.csv")
bulk_meta_tab %>% head()

# merge snRNA meta data with proteomics meta data -------------------------
meta_tab <- merge(snRNA_meta_tab[, c("Subject.ID", "Parent.Tumor.Segment.ID", "Current.Vial.Label")],
                  bulk_meta_tab[bulk_meta_tab$Type == "Tumor", c("Case.ID", "Specimen.Label")], 
                  by.x = c("Subject.ID"), by.y = c("Case.ID"), all.x = T)
meta_tab %>% head()
meta_tab <- meta_tab %>%
  rename(Case.ID = Subject.ID) %>%
  rename(Specimen.ID.bulk = Specimen.Label) %>%
  mutate(Specimen.ID.snRNA = str_replace_all(meta_tab$Current.Vial.Label, fixed(" "), "")) %>%
  select(-c("Current.Vial.Label"))


# Write meta data table ---------------------------------------------------
version_tmp <- 1
write.table(x = meta_tab, file = paste0(makeOutDir(), "meta_data", ".", format(Sys.Date(), "%Y%m%d") , ".v", version_tmp ,".tsv"), quote = F, sep = "\t", row.names = F)
