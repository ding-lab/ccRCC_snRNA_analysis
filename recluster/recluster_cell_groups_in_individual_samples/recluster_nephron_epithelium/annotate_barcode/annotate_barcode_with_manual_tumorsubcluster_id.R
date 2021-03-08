# Yige Wu @WashU March 2020
## annotate each barcode to the manual tumor subcluster number
## for each individual sample, input the tumor-cell-reclustered object and extract meta data
## then annotate with manual tumor subcluster number based on the seurat assigned tumor subcluster number

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
## input seurat object paths
srat_paths <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/recluster_nephron_epithelium_cells_in_individual_samples/20200225.v1/Seurat_Object_Paths.Malignant_Nephron_Epithelium20200225.v1.tsv", data.table = F)
## input the cell to cell type table
sratcluster2manualcluster_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/Tumor_Subcluster/ccRCC_snRNA_Downstream_Processing - Individual.TumorCluster2Cell_Type.20200223.v1.tsv", data.table = F)

# process by each aliquot ----------------------------------------------------
barcode_metadata_df <- NULL
for (aliquot_tmp in srat_paths$Aliquot) {
  ## input srat object
  srat_path <- srat_paths$Path_seurat_object[srat_paths$Aliquot == aliquot_tmp]
  srat <- readRDS(file = srat_path)
  
  ## extract current meta data
  barcode_metadata_tmp <- srat@meta.data
  barcode_metadata_tmp$barcode <- rownames(barcode_metadata_tmp)
  
  ## bind with the super table
  barcode_metadata_df <- rbind(barcode_metadata_tmp, barcode_metadata_df)
}
## add manual cluster
barcode2manualcluster_df <- merge(x = barcode_metadata_df, 
                             y = sratcluster2manualcluster_df %>%
                               select(Aliquot, Cluster, Cluster_Manual), 
                             by.x = c("orig.ident", "seurat_clusters"), 
                             by.y = c("Aliquot", "Cluster"), all.x = T)

## reformat
barcode2manualcluster_df <- barcode2manualcluster_df %>%
  select(orig.ident, barcode, seurat_clusters, Cluster_Manual) %>%
  rename(seurat_cluster_id = seurat_clusters) %>%
  rename(manual_cluster_id = Cluster_Manual)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "barcode2tumorsubclusterid.", run_id, ".tsv")
write.table(x = barcode2manualcluster_df, file = file2write, sep = '\t', quote = F, row.names = F)

  
  