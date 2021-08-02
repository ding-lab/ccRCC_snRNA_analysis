# Yige Wu @WashU Aug 2020
## 2020-10-27 new cell type correction, the cell type1 is no longer reliable
## 2020-11-19 corrected a bunch of fibroblasts and myofibroblasts cell type for snATAC datasets
## 2020-11-21 corrected fibroblasts and myofibroblasts cell type for all samples, switched to map the cell type detailed
## 2020-11-30 added a new group with detailed epithelial cell types

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the barcode-cell-type table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_corrected_by_individual_sample_inspection/20201121.v1/31Aliquot.Barcode2CellType.20201121.v1.tsv", data.table = F)
## input the newly added tumor sample
barcode2celltype_patch_df1 <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/individual_sample/map_celltype_for_C3N-00317-T1/20210305.v1/C3N-00317-T1.Barcode2CellType.20210305.v1.tsv")
barcode2celltype_patch_df2 <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/individual_sample/map_celltype_for_C3N-00437-T1/20210423.v1/C3N-00317-T1.Barcode2CellType.20210423.v1.tsv")
barcode2celltype_patch_df3 <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/individual_sample/map_celltype_for_C3N-00242-N/20210802.v1/C3N-00242-N.Barcode2CellType.20210802.v1.tsv")
barcode2celltype_patch_df4 <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/individual_sample/map_celltype_for_C3L-00079-N/20210802.v1/C3L-00079-N.Barcode2CellType.20210802.v1.tsv")

# group detailed immune cell types into major immune cell groups ----------
table(barcode2celltype_df$Cell_type.shorter)
barcode2celltype_df$Cell_group13 <- barcode2celltype_df$Cell_type.shorter
barcode2celltype_df$Cell_group13[barcode2celltype_df$Cell_group13 %in% c("Macrophages", "Macrophages proliferating", "TRM", "Monocytes")] <- "Macrophages"
barcode2celltype_df$Cell_group13[barcode2celltype_df$Cell_group13 %in% c("B-cells", "Plasma")] <- "B-cells"
barcode2celltype_df$Cell_group13[barcode2celltype_df$Cell_group13 %in% c("CD4 CTL", "CD4 T-cells", "CD4 T-cells activated", "CD4 T-cells naive", "Tregs")] <- "CD4+ T-cells"
barcode2celltype_df$Cell_group13[barcode2celltype_df$Cell_group13 %in% c("CD8 CTL", "CD8 CTL exhausted", "CD8 T-cells preexhausted")] <- "CD8+ T-cells"
barcode2celltype_df$Cell_group13[barcode2celltype_df$Cell_group13 %in% c("cDC", "pDC")] <- "DC"
barcode2celltype_df$Cell_group13[barcode2celltype_df$Cell_group13 %in% c("NK cells strong", "NK cells weak")] <- "NK cells"
barcode2celltype_df$Cell_group13[barcode2celltype_df$Cell_group13 %in% c("Basophils", "CD4/CD8 proliferating", "Mixed myeloid/lymphoid", "Mast cells")] <- "Immune others"
barcode2celltype_df$Cell_group13[barcode2celltype_df$Cell_group13 %in% c("Transitional cells", "Tumor-like cells", "EMT tumor cells")] <- "Tumor cells"
barcode2celltype_df$Cell_group13[barcode2celltype_df$Cell_group13 %in% c("Normal-like cells")] <- "Normal epithelial cells"

# make a new group for the transitional cells -----------------------------
barcode2celltype_df <- barcode2celltype_df %>%
  mutate(Cell_group14_w_transitional = ifelse(Cell_type.shorter == "EMT tumor cells", "EMT tumor cells", Cell_group13))

# make a new group for the epithelial cell types -----------------------------
barcode2celltype_df <- barcode2celltype_df %>%
  mutate(Cell_group_w_epithelialcelltypes = ifelse(Cell_type.shorter == "Normal epithelial cells", Cell_type.detailed, Cell_group13))

# rename other cell groups ------------------------------------------------
barcode2celltype_df <- barcode2celltype_df %>%
  dplyr::rename(Cell_group4 = Cell_group.shorter) %>%
  dplyr::rename(Cell_group5 = Cell_group.detailed)

# add patch ---------------------------------------------------------------
barcode2celltype_df <- rbind(barcode2celltype_patch_df2, barcode2celltype_patch_df1, barcode2celltype_df, barcode2celltype_patch_df3, barcode2celltype_patch_df4)
table(barcode2celltype_df$Cell_group13)
table(barcode2celltype_df$Cell_group14_w_transitional)
table(barcode2celltype_df$Cell_group_w_epithelialcelltypes)
table(barcode2celltype_df$Cell_group5)
table(barcode2celltype_df$Cell_group4)
table(barcode2celltype_df$orig.ident) %>% length()

# write output ------------------------------------------------------------
## final check up
nrow(barcode2celltype_df) # 161312
file2write <- paste0(dir_out, "35Aliquot.Barcode2CellType.", run_id, ".tsv")
write.table(x = barcode2celltype_df, file = file2write, quote = F, sep = "\t", row.names = F)





