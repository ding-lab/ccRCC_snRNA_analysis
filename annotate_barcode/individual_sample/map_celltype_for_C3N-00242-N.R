# Yige Wu @WashU Aug 2021
## make barcode to cell type mapping table for the integrated dataset
## just for normal epithelial cells

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

# input dependencies --------------------------------------------------
## input cluster2celltype
cluster2celltype_df <- readxl::read_excel(path = "./Resources/snRNA_Processed_Data/Cell_Type_Assignment/Individual_AllClusters/C3N-00242-N.20210802.xlsx")
## inpus seurat object
srat <- readRDS(file = "./Resources/snRNA_Processed_Data/scRNA_auto/outputs/CPT0014470002/pf2500_fmin200_fmax10000_cmin2500_cmax10000_mito_max0.1/CPT0014470002_processed.rds")
## input the barcod2celltype info to replace
barcode2celltype_replace_df <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/individual_sample/map_celltype_for_C3N-00242-N_immune_reclustered/20210802.v1/C3N-00242-N.Barcode2CellType.20210802.v1.tsv")
## input doublet prediction
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20210729.v1/scrublet.united_outputs.20210729.v1.tsv", data.table = F)

# subset ------------------------------------------------------------------
## subset to non-doublet
aliquot_show <- "C3N-00242-N"
scrublets_df <- barcode2scrublet_df %>%
  filter(Aliquot_WU == aliquot_show) %>%
  filter(predicted_doublet)
barcodes_keep <- rownames(srat@meta.data)
barcodes_keep <- barcodes_keep[!(barcodes_keep %in% scrublets_df$Barcodes)]
srat <- subset(x = srat, cells = barcodes_keep)

# map cell type to barcode ------------------------------------------------
barcode2cluster_df <- srat@meta.data 
#nrow(barcode2cluster_df)
barcode2cluster_df$individual_barcode <- rownames(barcode2cluster_df)
#unique(barcode2cluster_df$seurat_clusters)
barcode2cluster_df$seurat_clusters <- as.numeric(as.vector(barcode2cluster_df$seurat_clusters))
#unique(barcode2cluster_df$seurat_clusters)
barcode2celltype_df <- merge(barcode2cluster_df, 
                             cluster2celltype_df, 
                             by.x = c("orig.ident", "seurat_clusters"), by.y = c("Aliquot", "Cluster"), all.x = T)
## format
barcode2celltype_df <- barcode2celltype_df %>%
  mutate(Cell_type1 = NA) %>%
  mutate(Cell_type2 = NA) %>%
  mutate(Cell_type3 = NA) %>%
  mutate(Cell_type4 = NA) %>%
  mutate(Id_TumorManualCluster = NA) %>%
  mutate(Id_SeuratCluster = NA) %>%
  mutate(integrated_barcode = NA) %>%
  mutate(Comment = NA) %>%
  mutate(Cell_group13 = NA) %>%
  mutate(Cell_group14_w_transitional = NA) %>%
  mutate(Cell_group_w_epithelialcelltypes = NA) %>%
  select(orig.ident, Cell_type.shorter, Cell_type.detailed, 
         Cell_group4, Cell_group5, Cell_type1, Cell_type2, Cell_type3, Cell_type4,
         Id_TumorManualCluster, Id_SeuratCluster,
         individual_barcode, integrated_barcode, Comment,
         Cell_group13, Cell_group14_w_transitional, Cell_group_w_epithelialcelltypes)

# replace -----------------------------------------------------------------
## replace
barcode2celltype_df <- rbind(barcode2celltype_df %>%
                               filter(!(individual_barcode %in% barcode2celltype_replace_df$individual_barcode)),
                             barcode2celltype_replace_df)
table(barcode2celltype_df$Cell_group13)
table(barcode2celltype_df$Cell_type.shorter)
table(barcode2celltype_df$Cell_type.detailed)

# group detailed immune cell types into major immune cell groups ----------
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
table(barcode2celltype_df$Cell_group13)

# make a new group for the transitional cells -----------------------------
barcode2celltype_df <- barcode2celltype_df %>%
  mutate(Cell_group14_w_transitional = ifelse(Cell_type.shorter == "EMT tumor cells", "EMT tumor cells", Cell_group13))
table(barcode2celltype_df$Cell_group14_w_transitional)

# make a new group for the epithelial cell types -----------------------------
barcode2celltype_df <- barcode2celltype_df %>%
  mutate(Cell_group_w_epithelialcelltypes = ifelse(Cell_type.shorter == "Normal epithelial cells", Cell_type.detailed, Cell_group13))
table(barcode2celltype_df$Cell_group_w_epithelialcelltypes)


# write output ------------------------------------------------------------
file2write <- paste0(dir_out, aliquot_show, ".", "Barcode2CellType.", run_id, ".tsv")
write.table(x = barcode2celltype_df, file = file2write, quote = F, sep = "\t", row.names = F)
  