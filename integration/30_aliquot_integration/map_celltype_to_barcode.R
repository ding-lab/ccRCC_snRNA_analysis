# Yige Wu @WashU Feb 2020
## make barcode to cell type mapping table for the integrated dataset

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
## input barcode to cluster mapping table from 
all_integrated_barcode2cluster_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/fetch_data/20200212.v3/30_aliquot_integration.20200212.v3.umap_data.tsv", data.table = F)
## input cluster to cell type mapping table for all clusters in the integrated dataset
all_integrated_cluster2celltype_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/integration.allcluster2celltype.20200213.v3.tsv")
## input barcode to cluster mapping table for all immune cell integrated
immune_integrated_barcode2cluster_df <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_integrated_data/fetch_data/20200214.v1/integration.202002012.v3.immune_reclustered.20200213.v2.umap_data.20200214.v1.tsv", data.table = F)
## input cluster to cell type mapping table for immune clusters in the integrated dataset
immune_integrated_cluster2celltype_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/integration.immune.cluster2celltype.20200224.v1.csv", data.table = F)

# map cell types for immune cell group ----------------------------------------------------------
## merge barcode to cluster info with cluster to cell type info
immune_integrated_barcode2celltype_df <- merge(immune_integrated_barcode2cluster_df, immune_integrated_cluster2celltype_df,
                                                  by.x = c("ident"), by.y = c("Cluster"), all.x = T)
## format the barcode column and others to bind with other info later
immune_integrated_barcode2celltype_df <- immune_integrated_barcode2celltype_df %>%
  mutate(individual_barcode = str_split_fixed(string = barcode, pattern = "_", n = 2)[,1]) %>%
  mutate(Is_Normal_Nephron_Epithelium = F) %>%
  select(orig.ident, individual_barcode, Most_Enriched_Cell_Group, Most_Enriched_Cell_Type1, Most_Enriched_Cell_Type1, Most_Enriched_Cell_Type2, Most_Enriched_Cell_Type3, Most_Enriched_Cell_Type4, Is_Normal_Nephron_Epithelium)

# map cell types for cell groups other than immune----------------------------------------------------------
## merge barcode to cluster info with cluster to cell type info
all_integrated_barcode2celltype_df <- merge(all_integrated_barcode2cluster_df, all_integrated_cluster2celltype_df,
                                              by.x = c("ident"), by.y = c("Cluster"), all.x = T)
## format the barcode column and others to bind with other info later
all_integrated_barcode2celltype_df <- all_integrated_barcode2celltype_df %>%
  mutate(individual_barcode = str_split_fixed(string = barcode, pattern = "_", n = 2)[,1]) %>%
  mutate(Is_Normal_Nephron_Epithelium = ((ident == 17) & (Most_Enriched_Cell_Group == "Nephron_Epithelium"))) %>%
  select(orig.ident, individual_barcode, Most_Enriched_Cell_Group, Most_Enriched_Cell_Type1, Most_Enriched_Cell_Type1, Most_Enriched_Cell_Type2, Most_Enriched_Cell_Type3, Most_Enriched_Cell_Type4, Is_Normal_Nephron_Epithelium)

# merge the immune and non-immune info ------------------------------------
barcode2celltype_df <- rbind(immune_integrated_barcode2celltype_df,
                             all_integrated_barcode2celltype_df[all_integrated_barcode2celltype_df$Most_Enriched_Cell_Group != "Immune",])

# write output ------------------------------------------------------------
write.table(x = barcode2celltype_df, file = paste0(dir_out, "30_aliquot_integration.barcode2celltype.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)
