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
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input barcode to cluster mapping table from 
all_integrated_barcode2cluster_df <- fread(input = "./Resources/Analysis_Results/integration/30_aliquot_integration/fetch_data/20200212.v3/30_aliquot_integration.20200212.v3.umap_data.tsv", data.table = F)
## input cluster to cell type mapping table for all clusters in the integrated dataset
all_integrated_cluster2celltype_df <- fread(input = "./Resources/snRNA_Processed_Data/Cell_Type_Assignment/Integration_AllClusters/integration.allcluster2celltype.20200213.v3.tsv")
## input individual meta data
normal_reclustered_metadata_df <- fread(input = "./Resources/Analysis_Results/recluster/recluster_cell_groups_in_integrated_data/fetch_data/fetch_data_for_normal_epithelial_cells_local/20200406.v1/normal_epithelial_reclustered..metadata.20200406.v1.tsv", data.table = F)
## input individual cluster-to-cell type data
normal_reclustered_cluster2celltype_df <- fread(input = "./Resources/snRNA_Processed_Data/Cell_Type_Assignment/Normal_Epithelial_Cells/normal_epithelial_reclustered.cluster2celltype.20200406.v1.csv", data.table = F)

# label the normal epithelial cells in the integrated data ---------------------------------------
## merge barcode to cluster info with cluster to cell type info
all_integrated_barcode2celltype_df <- merge(all_integrated_barcode2cluster_df, all_integrated_cluster2celltype_df,
                                            by.x = c("ident"), by.y = c("Cluster"), all.x = T)
## format the barcode column and others to bind with other info later
all_integrated_barcode2celltype_df <- all_integrated_barcode2celltype_df %>%
  mutate(individual_barcode = str_split_fixed(string = barcode, pattern = "_", n = 2)[,1]) %>%
  mutate(Is_Normal_Nephron_Epithelium = ((ident == 17) & (Most_Enriched_Cell_Group == "Nephron_Epithelium"))) %>%
  mutate(integrated_barcode = barcode)
## get normal epithelial cells
normal_epithelial_barcode2celltype_df <- all_integrated_barcode2celltype_df %>%
  filter(Is_Normal_Nephron_Epithelium == T)

# map cell type based on normal-epithelial-reclustered cluster-based cell type -------------------------------------------------------
## add by cluster-based cell type to all barcodes
normal_reclustered_barcode2celltype_df <- merge(normal_reclustered_metadata_df, 
                                                normal_reclustered_cluster2celltype_df,
                                                by.x = c("seurat_clusters"), 
                                                by.y = c("Cluster"), all.x = T)

## merge with normal epithelial cell barcodes
normal_epithelial_barcode2celltype_df <- merge(normal_epithelial_barcode2celltype_df %>%
                                                 select(orig.ident, individual_barcode, integrated_barcode, Is_Normal_Nephron_Epithelium),
                                               normal_reclustered_barcode2celltype_df %>%
                                                 mutate(individual_barcode = str_split_fixed(string = normal_integrated_barcode, pattern = "_", n = 2)[,1]) %>%
                                                 select(orig.ident, individual_barcode, 
                                                        Most_Enriched_Cell_Group, Most_Enriched_Cell_Type1, Most_Enriched_Cell_Type2, Most_Enriched_Cell_Type3, Most_Enriched_Cell_Type4),
                                               by = c("orig.ident", "individual_barcode"),
                                               all.x = T)

## make detailed cell type
normal_epithelial_barcode2celltype_df$Cell_type.detailed <- normal_epithelial_barcode2celltype_df$Most_Enriched_Cell_Type2
normal_epithelial_barcode2celltype_df$Cell_type.detailed[normal_epithelial_barcode2celltype_df$Cell_type.detailed == ""] <- normal_epithelial_barcode2celltype_df$Most_Enriched_Cell_Type1[normal_epithelial_barcode2celltype_df$Cell_type.detailed == ""]
normal_epithelial_barcode2celltype_df$Cell_type.detailed[normal_epithelial_barcode2celltype_df$Cell_type.detailed == ""] <- normal_epithelial_barcode2celltype_df$Most_Enriched_Cell_Group[normal_epithelial_barcode2celltype_df$Cell_type.detailed == ""]
normal_epithelial_barcode2celltype_df$Cell_type.detailed[normal_epithelial_barcode2celltype_df$Cell_type.detailed == "Thick ascending limb"] <- "Loop of Henle"
## make detailed cell type
normal_epithelial_barcode2celltype_df$Cell_type.shorter <- "Normal epithelial cells"

## add columns
normal_epithelial_barcode2celltype_df <- normal_epithelial_barcode2celltype_df %>%
  mutate(is_malignant = F) %>%
  mutate(manual_tumorsubcluster_id = NA) %>%
  select(orig.ident, individual_barcode, integrated_barcode, 
         Most_Enriched_Cell_Group, is_malignant,
         Cell_type.shorter, Cell_type.detailed, 
         Most_Enriched_Cell_Type1, Most_Enriched_Cell_Type2, Most_Enriched_Cell_Type3, Most_Enriched_Cell_Type4, manual_tumorsubcluster_id, Is_Normal_Nephron_Epithelium)

# write output ------------------------------------------------------------
write.table(x = normal_epithelial_barcode2celltype_df, file = paste0(dir_out, "normal_epithelial_reclustered.barcode2celltype.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)
