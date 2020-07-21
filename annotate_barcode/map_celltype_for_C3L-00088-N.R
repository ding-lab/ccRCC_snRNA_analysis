# Yige Wu @WashU Apr 2020
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
cluster2celltype_df <- fread(data.table = F, input = "./Resources/snRNA_Processed_Data/Cell_Type_Assignment/Individual_AllClusters/ccRCC_snRNA_Downstream_Processing - Individual.AllCluster2Cell_Type.20200630.v1.tsv")
## inpus seurat object
srat <- readRDS(file = "./Resources/snRNA_Processed_Data/scRNA_auto/outputs/CPT0000890002/pf1000_fmin200_fmax10000_cmin1000_cmax10000_mito_max0.1/CPT0000890002_processed.rds")

# map cell type to barcode ------------------------------------------------
barcode2cluster_df <- srat@meta.data 
barcode2cluster_df <- as.data.frame(barcode2cluster_df)
barcode2cluster_df$individual_barcode <- rownames(barcode2cluster_df)
barcode2cluster_df$seurat_clusters <- as.numeric(barcode2cluster_df$seurat_clusters)
barcode2celltype_df <- merge(barcode2cluster_df, 
                             cluster2celltype_df %>%
                               filter(Aliquot == "CPT0000890002"), 
                             by.x = c("orig.ident", "seurat_clusters"), by.y = c("Aliquot", "Cluster"))
## add cell type detailed
barcode2celltype_df$Cell_type.detailed <- barcode2celltype_df$Most_Enriched_Cell_Type2
barcode2celltype_df$Cell_type.detailed[barcode2celltype_df$Cell_type.detailed == ""] <- barcode2celltype_df$Most_Enriched_Cell_Type1[barcode2celltype_df$Cell_type.detailed == ""]
barcode2celltype_df$Cell_type.detailed %>% unique()
barcode2celltype_df$Cell_type.detailed[barcode2celltype_df$Cell_type.detailed == "Thick ascending limb"] <- "Loop of Henle"

## add cell type shorter
barcode2celltype_df$Cell_type.shorter <- barcode2celltype_df$Cell_type.detailed
barcode2celltype_df$Cell_type.shorter[barcode2celltype_df$Most_Enriched_Cell_Group == "Nephron_Epithelium"] <- "Normal epithelial cells"
## format
barcode2celltype_df <- barcode2celltype_df %>%
  mutate(integrated_barcode = NA) %>%
  mutate(Id_TumorManualCluster = NA) %>%
  mutate(Cell_group = ifelse(Most_Enriched_Cell_Group == "Nephron_Epithelium", "Normal epithelial cells", Most_Enriched_Cell_Group)) %>%
  select(orig.ident, individual_barcode, integrated_barcode, Most_Enriched_Cell_Group, Cell_type.shorter, Cell_type.detailed, Most_Enriched_Cell_Type1, Most_Enriched_Cell_Type2, Most_Enriched_Cell_Type3, Most_Enriched_Cell_Type4, Id_TumorManualCluster, Cell_group)
  
# write output ------------------------------------------------------------
nrow(barcode2celltype_df)
file2write <- paste0(dir_out, "Barcode2CellType.", "C3L-00088-N.", run_id, ".tsv")
write.table(x = barcode2celltype_df, file = file2write, quote = F, sep = "\t", row.names = F)
  