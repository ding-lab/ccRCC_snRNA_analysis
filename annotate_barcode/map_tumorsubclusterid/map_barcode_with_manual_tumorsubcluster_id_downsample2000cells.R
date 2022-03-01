# Yige Wu @WashU Dec 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "Seurat"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- "cutoff50cells"
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory to be in the same structure as the code
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input barcode2seurat cluster info
barcode2seuratcluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/recluster/recluster_tumorcells/downsample_fixednumber_and_recluster_tumor_cells_in_selected_samples_katmai/20220222.v1/UMAPData.2000TumorCellReclustered.20220222.v1.tsv")
## input the cell to cell type table
sratcluster2manualcluster_df <- readxl::read_xlsx(path = "./Resources/snRNA_Processed_Data/Tumor_Subclusters/Individual.Downsample2000cells.TumorSeuratCluster2Manual.20220223.xlsx")
## input meta data
### accidently use the wrong sample id
metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210809.v1/meta_data.20210809.v1.tsv")

# get barcode2manualcluster -----------------------------------------------
barcode2seuratcluster_df$easy_id <- mapvalues(x = barcode2seuratcluster_df$orig.ident, from = metadata_df$Aliquot.snRNA, to = as.vector(metadata_df$Aliquot.snRNA.WU))
barcode2manualcluster_df <- merge(x = barcode2seuratcluster_df %>%
                                    filter(easy_id != "C3L-00359-T1") %>%
                                    rename(id_seurat_cluster = seurat_clusters), 
                                  y = sratcluster2manualcluster_df %>%
                                    mutate(id_manual_cluster_w0 = ifelse(id_manual_cluster_w0 == "NA", NA, id_manual_cluster_w0)) %>%
                                    mutate(id_manual_cluster_w0 = as.numeric(id_manual_cluster_w0)) %>%
                                    select(Aliquot, id_seurat_cluster, id_manual_cluster_w0, Comment), 
                                  by.x = c("orig.ident", "id_seurat_cluster"), 
                                  by.y = c("Aliquot", "id_seurat_cluster"), all.x = T)
barcode2manualcluster_df <- barcode2manualcluster_df %>%
  mutate(Cluster_Name = paste0(easy_id, "_C", id_manual_cluster_w0+1))

cellnumber_percluster_df <- barcode2manualcluster_df %>%
  # filter(!is.na(id_manual_cluster_w0)) %>%
  group_by(Cluster_Name) %>%
  summarise(Freq = n())

barcode2manualcluster_df <- barcode2manualcluster_df %>%
  mutate(Cluster_Name.cutoff50cells = ifelse(Cluster_Name %in% cellnumber_percluster_df$Cluster_Name[cellnumber_percluster_df$Freq >= 50], Cluster_Name, paste0(easy_id, "_CNA")))

barcode2manualcluster_df %>%
  select(Cluster_Name, id_manual_cluster_w0) %>%
  unique()

unique(barcode2manualcluster_df$Cluster_Name[!is.na(barcode2manualcluster_df$id_manual_cluster_w0)])
## 98

unique(barcode2manualcluster_df$Cluster_Name.cutoff50cells[!grepl(pattern = "CNA", x = barcode2manualcluster_df$Cluster_Name.cutoff50cells)])
## 83

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Barcode2TumorSubclusterId.", run_id, ".tsv")
write.table(x = barcode2manualcluster_df, file = file2write, sep = '\t', quote = F, row.names = F)
