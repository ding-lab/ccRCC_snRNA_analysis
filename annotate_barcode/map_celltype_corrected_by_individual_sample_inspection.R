# Yige Wu @WashU Aug 2020
## make barcode to cell type mapping table for cell types changed based on individual sample inspection

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(dplyr)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the integrated-data-based barcode to cell type table
barcode2celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_to_all_cells_with_patch/20200904.v1/31AliquotIntegration.Barcode2CellType.TumorManualCluster.20200904.v1.tsv")
## input barcode to individual cluster id 
barcode2metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/data_summary/fetch_data/fetch_data_by_individual_sample/20200717.v1/Barcode2MetaData.20200717.v1.tsv")
## input corrected cell type
celltypecorrected_df <- readxl::read_excel(path = "./Resources/snRNA_Processed_Data/Cell_Type_Assignment/Individual_AllClusters/Cells_BySampleByClusterByCellTypeShorter.Over50.20200917.xlsx", sheet = "Sheet1")

# merge info -------------------------------------
## merge cell type with seurat cluster
merged_df <- merge(x = barcode2celltype_df, barcode2metadata_df, 
                   by.x = c("orig.ident", "individual_barcode"), by.y = c("aliquot", "individual_barcode"), all.x = T)
# which(is.na(merged_df$seurat_cluster_id))
## change column names
merged_df <- merged_df %>%
  dplyr::rename(Cell_group.shorter = Most_Enriched_Cell_Group) %>%
  dplyr::rename(Cell_group.detailed = Cell_group) %>%
  dplyr::rename(Cell_type1 = Most_Enriched_Cell_Type1) %>%
  dplyr::rename(Cell_type2 = Most_Enriched_Cell_Type2) %>%
  dplyr::rename(Cell_type3 = Most_Enriched_Cell_Type3) %>%
  dplyr::rename(Cell_type4 = Most_Enriched_Cell_Type4) %>%
  dplyr::rename(Id_SeuratCluster = seurat_cluster_id) %>%
  dplyr::select(orig.ident,
         Cell_type.shorter, Cell_type.detailed, 
         Cell_group.shorter, Cell_group.detailed,
         Cell_type1, Cell_type2, Cell_type3, Cell_type4, 
         Id_TumorManualCluster, Id_SeuratCluster,
         individual_barcode, integrated_barcode) %>%
  dplyr::mutate(Id_Mapping_Cells = paste0(orig.ident, "_", Id_SeuratCluster, "_", Cell_type.shorter))
## add mapping id to the corrected cell type map
celltypecorrected_df <- celltypecorrected_df %>%
  dplyr::mutate(Id_Mapping_Cells = paste0(aliquot, "_", seurat_cluster_id, "_", Cell_type.shorter.original))
celltypecorrected_df$Cell_type1[is.na(celltypecorrected_df$Cell_type1)] <- ""
## correct cell type
### Cell_type.shorter
temp <- mapvalues(x = merged_df$Id_Mapping_Cells, from = celltypecorrected_df$Id_Mapping_Cells, to = as.vector(celltypecorrected_df$Cell_type.shorter))
idx_kepp <- (temp == merged_df$Id_Mapping_Cells)
temp[idx_kepp] <- merged_df$Cell_type.shorter[idx_kepp]
merged_df$Cell_type.shorter <- temp
table(merged_df$Cell_type.shorter)
### Cell_type.detailed
temp <- mapvalues(x = merged_df$Id_Mapping_Cells, from = celltypecorrected_df$Id_Mapping_Cells, to = as.vector(celltypecorrected_df$Cell_type.detailed))
idx_kepp <- (temp == merged_df$Id_Mapping_Cells)
temp[idx_kepp] <- merged_df$Cell_type.detailed[idx_kepp]
merged_df$Cell_type.detailed <- temp
table(merged_df$Cell_type.detailed)
### Cell_group.shorter
temp <- mapvalues(x = merged_df$Id_Mapping_Cells, from = celltypecorrected_df$Id_Mapping_Cells, to = as.vector(celltypecorrected_df$Cell_group.shorter))
idx_kepp <- (temp == merged_df$Id_Mapping_Cells)
temp[idx_kepp] <- merged_df$Cell_group.shorter[idx_kepp]
merged_df$Cell_group.shorter <- temp
table(merged_df$Cell_group.shorter)
### Cell_group.detailed
temp <- mapvalues(x = merged_df$Id_Mapping_Cells, from = celltypecorrected_df$Id_Mapping_Cells, to = as.vector(celltypecorrected_df$Cell_group.detailed))
idx_kepp <- (temp == merged_df$Id_Mapping_Cells)
temp[idx_kepp] <- merged_df$Cell_group.detailed[idx_kepp]
merged_df$Cell_group.detailed <- temp
table(merged_df$Cell_group.detailed)
### Cell_type1
temp <- mapvalues(x = merged_df$Id_Mapping_Cells, from = celltypecorrected_df$Id_Mapping_Cells, to = as.vector(celltypecorrected_df$Cell_type1))
idx_kepp <- (temp == merged_df$Id_Mapping_Cells)
temp[idx_kepp] <- merged_df$Cell_type1[idx_kepp]
merged_df$Cell_type1 <- temp
table(merged_df$Cell_type1)
## add comment
temp <- mapvalues(x = merged_df$Id_Mapping_Cells, from = celltypecorrected_df$Id_Mapping_Cells, to = as.vector(celltypecorrected_df$Comment))
idx_kepp <- (temp == merged_df$Id_Mapping_Cells)
temp[idx_kepp] <- ""
merged_df$Comment <- temp
table(merged_df$Comment)
## remove column
merged_df <- merged_df %>%
  dplyr::select(-Id_Mapping_Cells)
# write output ------------------------------------------------------------
## final check up
nrow(merged_df) # 138547
nrow(barcode2celltype_df) # 138547
## 94846 tumor cells, 989 unknown cells, 213 tumor-like epithelial cells
file2write <- paste0(dir_out, "31Aliquot.Barcode2CellType.", run_id, ".tsv")
write.table(x = merged_df, file = file2write, quote = F, sep = "\t", row.names = F)

