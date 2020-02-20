# Yige Wu @WashU Feb 2020
## for calculating the tumor content, by taking out the normal-like barcodes
## the rest mostly show copy number alterations by manual inspection

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
## input cluster to cell type table (most updated)
cluster2celltype_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - Individual.AllCluster2Cell_Type.20200207.v2.tsv", data.table = F)
## input barcodes co-clustered with NAT in the integration
c17_barcodes_clean <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/other_plotting/umap_normal_like_cells_by_aliquot/20200218.v1/30_aliquot_integration.20200212.v3.Normal_Like_Cluster17_Barcodes.20200218.v1.tsv", data.table = F)
## input meta data file to distinguish origin tumor piece and additional tumor piece
meta_data_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/make_meta_data/20191105.v1/meta_data.20191105.v1.tsv", data.table = F)

# sum the nephron epithelium and take out normal-like cells ---------------
tumor_barcode_df <- NULL
for (snRNA_aliquot_id_tmp in unique(cluster2celltype_df$Aliquot)) {
  ## input the barcode to cluster table from infercnv input files
  path_file <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/InferCNV/inputs/annotations_file/Individual.20200207.v1/", snRNA_aliquot_id_tmp, ".Barcode_Annotation.txt")
  barcode_cluster_tmp <- fread(input = path_file, data.table = F)
  
  ## count the total barcodes in this aliquot
  total_barcodes <- nrow(barcode_cluster_tmp)
  
  ## filter down to only barcodes assigned as nephron epithelium
  nephron_epithelium_clusters <- cluster2celltype_df$Cluster[cluster2celltype_df$Aliquot == snRNA_aliquot_id_tmp & cluster2celltype_df$Most_Enriched_Cell_Group == "Nephron_Epithelium"]
  barcode_cluster_tmp <- barcode_cluster_tmp %>%
    rename(cluster = V2) %>%
    rename(barcode = V1) %>%
    filter(cluster %in% nephron_epithelium_clusters)
  
  ## take out barcodes clustered with NAT from the integration
  normal_like_barcodes <- c17_barcodes_clean$barcode_clean[c17_barcodes_clean$aliquot == snRNA_aliquot_id_tmp]
  if (length(normal_like_barcodes) > 0) {
    barcode_cluster_tmp <- barcode_cluster_tmp %>%
      filter(!(barcode %in% normal_like_barcodes))
  }
  
  ## count barcodes for remaining tumor cells
  tumor_barcodes <- nrow(barcode_cluster_tmp)
  
  ## assemble the two barcode number into a data frame
  tumor_barcode_df_tmp <- data.frame(Aliquot = snRNA_aliquot_id_tmp,
                                     Num_Tumor_Barcodes = tumor_barcodes,
                                     Num_Total_Barcodes = total_barcodes)
  tumor_barcode_df <- rbind(tumor_barcode_df_tmp, tumor_barcode_df)
}
## calculate tumor fractin
tumor_barcode_df <- tumor_barcode_df %>%
  mutate(Frac_Tumor_Barcodes = Num_Tumor_Barcodes/Num_Total_Barcodes)
## merge with meta data to distinguish original tumor piece and additional ones
tumor_barcode_df <- merge(tumor_barcode_df, meta_data_df, by.x = c("Aliquot"), by.y = c("Aliquot.snRNA"), all.x = T)

# write.table -------------------------------------------------------------
write.table(x = tumor_barcode_df, file = paste0(dir_out, "Tumor_Cell_Fraction.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)

