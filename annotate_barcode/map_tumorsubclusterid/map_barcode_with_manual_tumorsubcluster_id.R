# Yige Wu @WashU Jun 2020
## annotate each barcode to the manual tumor subcluster number (raw)

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
## input barcode2seurat cluster info
barcode2seuratcluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/fetch_data/fetchdata_individual_tumorcellreclustered_on_katmai/20201127.v1/MetaData_TumorCellOnlyReclustered.20201127.v1.tsv")
## input the cell to cell type table
sratcluster2manualcluster_df <- readxl::read_xlsx(path = "./Resources/snRNA_Processed_Data/Tumor_Subclusters/Individual.TumorSeuratCluster2Manual.20201127.v1.xlsx")

# get barcode2manualcluster -----------------------------------------------
barcode2manualcluster_df <- merge(x = barcode2seuratcluster_df %>%
                                    mutate(id_seurat_cluster = seurat_clusters) %>%
                                    mutate(barcode = barcode_tumorcellreclustered) %>%
                                    select(easy_id, orig.ident, id_seurat_cluster, barcode), 
                                  y = sratcluster2manualcluster_df %>%
                                    select(Aliquot, id_seurat_cluster, id_manual_cluster_w0, Comment), 
                                  by.x = c("orig.ident", "id_seurat_cluster"), 
                                  by.y = c("Aliquot", "id_seurat_cluster"), all.x = T)
barcode2manualcluster_df <- barcode2manualcluster_df %>%
  mutate(Cluster_Name = paste0(easy_id, "_C", id_manual_cluster_w0+1))

length(unique(barcode2manualcluster_df$Cluster_Name[!is.na(barcode2manualcluster_df$id_manual_cluster_w0)]))
## 134
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Barcode2TumorSubclusterId.", run_id, ".tsv")
write.table(x = barcode2manualcluster_df, file = file2write, sep = '\t', quote = F, row.names = F)
