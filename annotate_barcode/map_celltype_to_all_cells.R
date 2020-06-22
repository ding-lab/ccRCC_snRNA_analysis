# Yige Wu @WashU Apr 2020
## make barcode to cell type mapping table for the integrated dataset

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
## input Alla's immune cell type assignment
immune_barcode2celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_to_immune_cells/20200616.v1/Barcode2ImmuneCellType.20200616.v1.tsv")
## input tumor subclustering cell type assignment: 93277 cells
tumor_barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_barcode_with_manual_tumorsubcluster_id/20200616.v1/Barcode2TumorSubclusterId.20200616.v1.tsv", data.table = F)
## input barcode-to-cell-type table
normal_epithelial_barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_to_normal_epithelial_cells/20200410.v1/normal_epithelial_reclustered.barcode2celltype.20200410.v1.tsv", data.table = F)

# format normal epithelial cells ------------------------------------------
normal_epithelial_barcode2celltype_df <- normal_epithelial_barcode2celltype_df %>%
  mutate(Id_TumorManualCluster = NA) %>%
  select(orig.ident, individual_barcode, integrated_barcode, 
         Most_Enriched_Cell_Group,
         Cell_type.shorter, Cell_type.detailed, 
         Most_Enriched_Cell_Type1, Most_Enriched_Cell_Type2, Most_Enriched_Cell_Type3, Most_Enriched_Cell_Type4, Id_TumorManualCluster)

# map stromal cells -------------------------------------------------------
stroma_barcode2celltype_df <- all_integrated_barcode2celltype_df %>%
  filter(Most_Enriched_Cell_Group == "Stroma") %>%
  mutate(Cell_type.shorter = Most_Enriched_Cell_Type1) %>%
  mutate(Cell_type.detailed = Most_Enriched_Cell_Type1)
## see what needs to be corrected
table(stroma_barcode2celltype_df$Most_Enriched_Cell_Type1)
table(stroma_barcode2celltype_df$Most_Enriched_Cell_Type2)
## make short cell type names for myofibroblasts
stroma_barcode2celltype_df$Cell_type.shorter[stroma_barcode2celltype_df$Most_Enriched_Cell_Type2 == "Myofibroblasts"] <- "Myofibroblasts"
stroma_barcode2celltype_df$Cell_type.detailed[stroma_barcode2celltype_df$Most_Enriched_Cell_Type2 == "Myofibroblasts"] <- "Myofibroblasts"
## add columns
stroma_barcode2celltype_df <- stroma_barcode2celltype_df %>%
  mutate(Id_TumorManualCluster = NA) %>%
  select(orig.ident, individual_barcode, integrated_barcode, 
         Most_Enriched_Cell_Group,
         Cell_type.shorter, Cell_type.detailed, 
         Most_Enriched_Cell_Type1, Most_Enriched_Cell_Type2, Most_Enriched_Cell_Type3, Most_Enriched_Cell_Type4, Id_TumorManualCluster)

# map unknown cells -------------------------------------------------------
assigned_integrated_barcodes <- c(immune_integrated_barcode2celltype_df$integrated_barcode,
                                  tumor_barcode2celltype_df$integrated_barcode,
                                  stroma_barcode2celltype_df$integrated_barcode,
                                  normal_epithelial_barcode2celltype_df$integrated_barcode)
length(assigned_integrated_barcodes)
nrow(all_integrated_barcode2celltype_df)

unknown_barcode2celltype_df <- all_integrated_barcode2celltype_df %>%
  filter(!(integrated_barcode %in% assigned_integrated_barcodes)) %>%
  mutate(Id_TumorManualCluster = NA) %>%
  mutate(Cell_type.shorter = "Unknown") %>%
  mutate(Cell_type.detailed = "Unknown")
unknown_barcode2celltype_df$Most_Enriched_Cell_Group <- "Unknown"
unknown_barcode2celltype_df$Most_Enriched_Cell_Type1 <- ""
unknown_barcode2celltype_df$Most_Enriched_Cell_Type2 <- ""
unknown_barcode2celltype_df$Most_Enriched_Cell_Type3 <- ""
unknown_barcode2celltype_df$Most_Enriched_Cell_Type4 <- ""

unknown_barcode2celltype_df <- unknown_barcode2celltype_df %>%
  select(orig.ident, individual_barcode, integrated_barcode, 
         Most_Enriched_Cell_Group,
         Cell_type.shorter, Cell_type.detailed, 
         Most_Enriched_Cell_Type1, Most_Enriched_Cell_Type2, Most_Enriched_Cell_Type3, Most_Enriched_Cell_Type4, Id_TumorManualCluster)

# merge the immune and tumor cells and other info ------------------------------------
barcode2celltype_df <- rbind(immune_barcode2celltype_df %>%
                               mutate(Cell_group = ifelse(Most_Enriched_Cell_Group == "Immune", "Immune", "Unknown")),
                             tumor_barcode2celltype_df %>%
                               mutate(Cell_group = ifelse(Cell_type.shorter == "Tumor cells", "Tumor cells", "Unknown")),
                             stroma_barcode2celltype_df %>%
                               mutate(Cell_group = "Stroma"),
                             normal_epithelial_barcode2celltype_df %>%
                               mutate(Cell_group = ifelse(Most_Enriched_Cell_Group == "Nephron_Epithelium", "Normal epithelial cells", "Unknown")),
                             unknown_barcode2celltype_df %>%
                               mutate(Cell_group = "Unknown"))
nrow(barcode2celltype_df)
table(barcode2celltype_df$Cell_type.shorter)
# write output ------------------------------------------------------------
write.table(x = barcode2celltype_df, file = paste0(dir_out, "30AliquotIntegration.Barcode2CellType.TumorManualCluster.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)
