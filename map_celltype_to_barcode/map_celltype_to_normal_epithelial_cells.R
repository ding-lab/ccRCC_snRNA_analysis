# Yige Wu @WashU Apr 2020
## make barcode to cell type mapping table for the integrated dataset
## just for normal epithelial cells

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
## input Alla's immune cell type assignment
immune_integrated_barcode2celltype_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/Immune/30_aliquot_integration.barcode2celltype.20200330.AK.v1.tsv", data.table = F)
## input tumor subclustering cell type assignment
barcode2manualtumorcluster_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/annotate_barcode/annotate_barcode_with_manual_tumorsubcluster_id/20200324.v1/barcode2tumorsubclusterid.20200324.v1.tsv", data.table = F)

# label the normal epithelial cells ---------------------------------------
## merge barcode to cluster info with cluster to cell type info
all_integrated_barcode2celltype_df <- merge(all_integrated_barcode2cluster_df, all_integrated_cluster2celltype_df,
                                            by.x = c("ident"), by.y = c("Cluster"), all.x = T)
## format the barcode column and others to bind with other info later
all_integrated_barcode2celltype_df <- all_integrated_barcode2celltype_df %>%
  mutate(individual_barcode = str_split_fixed(string = barcode, pattern = "_", n = 2)[,1]) %>%
  mutate(Is_Normal_Nephron_Epithelium = ((ident == 17) & (Most_Enriched_Cell_Group == "Nephron_Epithelium"))) %>%
  mutate(integrated_barcode = barcode)

# map normal epithelial cells -------------------------------------------------------
normal_epithelial_barcode2celltype_df <- all_integrated_barcode2celltype_df %>%
  filter(Is_Normal_Nephron_Epithelium == T) %>%
  mutate(Cell_type.shorter = Most_Enriched_Cell_Type1) %>%
  mutate(Cell_type.detailed = Most_Enriched_Cell_Type2)
table(normal_epithelial_barcode2celltype_df$Most_Enriched_Cell_Type1)
table(normal_epithelial_barcode2celltype_df$Most_Enriched_Cell_Type2)
## add columns
normal_epithelial_barcode2celltype_df <- normal_epithelial_barcode2celltype_df %>%
  mutate(is_malignant = F) %>%
  mutate(manual_tumorsubcluster_id = NA) %>%
  select(orig.ident, individual_barcode, integrated_barcode, 
         Most_Enriched_Cell_Group, is_malignant,
         Cell_type.shorter, Cell_type.detailed, 
         Most_Enriched_Cell_Type1, Most_Enriched_Cell_Type2, Most_Enriched_Cell_Type3, Most_Enriched_Cell_Type4, manual_tumorsubcluster_id, Is_Normal_Nephron_Epithelium)

# write output ------------------------------------------------------------
write.table(x = barcode2celltype_df, file = paste0(dir_out, "30_aliquot_integration.barcode2celltype.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)
