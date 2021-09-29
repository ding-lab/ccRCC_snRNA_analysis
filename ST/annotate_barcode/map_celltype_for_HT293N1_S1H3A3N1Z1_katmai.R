# Yige Wu @WashU Sep 2020

# set up libraries and output directory -----------------------------------
## getting the path to the current script
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
path_this_script <- thisFile()
## set working directory
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input seurat object
srat <- readRDS(file = "../ccRCC_ST/Processed_Data/Seurat/Outputs/TWFU-HT293N1-S1H3A3N1Z1_1Bmn1_1_5.0/TWFU-HT293N1-S1H3A3N1Z1_1Bmn1_1_5.0_processed_multiomic.rds")
## input cluster2celltype
cluster2celltype_df <- readxl::read_excel(path = "./Resources/snRNA_Processed_Data/Cell_Type_Assignment/Individual_AllClusters/HT293N1_S1H3A3N1Z1.xlsx")

# subset ------------------------------------------------------------------
## subset to non-doublet
aliquot_show <- "HT293N1_S1H3A3N1Z1"

# map cell type to barcode ------------------------------------------------
barcode2cluster_df <- FetchData(object = srat, vars = c("UMAP_1", "UMAP_2", "orig.ident", "seurat_clusters"))
#nrow(barcode2cluster_df)
barcode2cluster_df$individual_barcode <- rownames(barcode2cluster_df)
#unique(barcode2cluster_df$seurat_clusters)
barcode2cluster_df$seurat_clusters <- as.numeric(as.vector(barcode2cluster_df$seurat_clusters))
#unique(barcode2cluster_df$seurat_clusters)
barcode2celltype_df <- merge(barcode2cluster_df, 
                             cluster2celltype_df, 
                             by.x = c("seurat_clusters"), by.y = c("Cluster"), all.x = T)
## format
barcode2celltype_df <- barcode2celltype_df %>%
  select(orig.ident, Cell_type.shorter, Cell_type.detailed, 
         Cell_group4, Cell_group5,
         individual_barcode, Comment)

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
## make a new group for the transitional cells
barcode2celltype_df <- barcode2celltype_df %>%
  mutate(Cell_group14_w_transitional = ifelse(Cell_type.shorter == "EMT tumor cells", "EMT tumor cells", Cell_group13))
table(barcode2celltype_df$Cell_group14_w_transitional)
## make a new group for the epithelial cell types
barcode2celltype_df <- barcode2celltype_df %>%
  mutate(Cell_group_w_epithelialcelltypes = ifelse(Cell_type.shorter == "Normal epithelial cells", Cell_type.detailed, Cell_group13))
table(barcode2celltype_df$Cell_group_w_epithelialcelltypes)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, aliquot_show, ".Barcode2CellType.", run_id, ".tsv")
write.table(x = barcode2celltype_df, file = file2write, quote = F, sep = "\t", row.names = F)
rownames(barcode2celltype_df) <- barcode2celltype_df$individual_barcode
srat@meta.data <- barcode2celltype_df
file2write <- paste0(dir_out, aliquot_show, ".multiomic.celltypeannotated.", run_id, ".RDS")
saveRDS(object = srat, file = file2write, compress = T)
