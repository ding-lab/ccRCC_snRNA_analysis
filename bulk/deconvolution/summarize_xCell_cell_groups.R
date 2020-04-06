# Yige Wu @WashU March 2020
## running on local
## for summarizing xCell-estimated fraction of different cell types from CPTAC ccRCC manuscript
## the input xcell result include 103 CCRCC Tumor Samples and 72 NAT Samples
## so a lot of NAT samples have 0 for immune cells

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
## input xcell results
xcell_value_wide_df <- readxl::read_excel("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC Manuscript/CPTAC3-ccRCC-SupplementaryTables_Final/Table S7.xlsx", sheet = "xCell Signatures", skip = 2)
# because there is a row in character format so the rest of the data frame is in character format
# need to change it to double format
xcell_value_wide_df <- xcell_value_wide_df %>%
  filter(Samples != "Immune Group")
# xcell_value_wide_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Bulk_Processed_Data/xCell/outputs/xCell_CCRCC_RNA-Seq_Expr_WashU_FPKM_UQ_xCell_1007022020.txt")
# xcell_value_wide_df <- xcell_value_wide_df %>%
#   rename(Samples = V1)
## input table of cell family from xCell
xcell_cell_family_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Bulk_Processed_Data/xCell/scripts/xCell/Dev_scripts/cells_families.txt", data.table = F)
## set group "patient" cells to extract
parent_cells2process <- c("Endothelial cells", "B-cells", "CD4+ T-cells", "CD8+ T-cells", "NK cells", "Fibroblasts", "Macrophages", "DC", "Epithelial cells")
## input id meta data with Aliquot ids for snRNA
id_metadata_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/make_meta_data/20191105.v1/meta_data.20191105.v1.tsv", data.table = F)

# transform xCell values to numeric format --------------------------------------------------

## transform wide data frame to long data frame 
xcell_value_long_df <- melt(xcell_value_wide_df, id.vars = c("Samples"))
colnames(xcell_value_long_df) <- c("Cells", "Aliquot.pro", "Value")
## make value numeric
xcell_value_long_df$Value <- as.numeric(as.vector(xcell_value_long_df$Value))
xcell_value_long_df %>% head()

# filter by cells column -----------------------------------
xcell_value_long_filtered_df <- xcell_value_long_df %>%
  filter(Cells %in% parent_cells2process) %>%
  rename(Cellgroup.xCell = Cells) %>%
  rename(Value.xCell = Value)
xcell_value_long_filtered_df$Aliquot.snRNA <- mapvalues(x = xcell_value_long_filtered_df$Aliquot.pro, from = id_metadata_df$Aliquot.bulk, to = as.vector(id_metadata_df$Aliquot.snRNA))
xcell_value_long_filtered_df$Aliquot.pro <- as.vector(xcell_value_long_filtered_df$Aliquot.pro)
xcell_value_long_filtered_df$Aliquot.snRNA[xcell_value_long_filtered_df$Aliquot.snRNA == xcell_value_long_filtered_df$Aliquot.pro] <- NA

# transform long data frame to wide data frame ----------------------------
xcell_value_wide_filtered_df <- dcast(data = xcell_value_long_filtered_df, formula = Aliquot.pro ~ Cellgroup.xCell, value.var = "Value.xCell")

# add snRNA aliquot ids ---------------------------------------------------
xcell_value_wide_filtered_df$Aliquot.snRNA <- mapvalues(x = xcell_value_wide_filtered_df$Aliquot.pro, from = id_metadata_df$Aliquot.bulk, to = as.vector(id_metadata_df$Aliquot.snRNA))
xcell_value_wide_filtered_df$Aliquot.pro <- as.vector(xcell_value_wide_filtered_df$Aliquot.pro)
xcell_value_wide_filtered_df$Aliquot.snRNA[xcell_value_wide_filtered_df$Aliquot.snRNA == xcell_value_wide_filtered_df$Aliquot.pro] <- NA

# write output -------------------------------------------------------------------
file2write <- paste0(dir_out, "xcell_parent_cells_fraction.long.", run_id, ".tsv")
write.table(x = xcell_value_long_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)

file2write <- paste0(dir_out, "xcell_parent_cells_fraction.wide.", run_id, ".tsv")
write.table(x = xcell_value_wide_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
