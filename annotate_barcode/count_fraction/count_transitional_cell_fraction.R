# Yige Wu @WashU Apr 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 3
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input cell type per barcode table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups/20200917.v2/31Aliquot.Barcode2CellType.20200917.v2.tsv", data.table = F)
## input seurat paths
paths_srat_df <- fread(data.table = F, input = "./Resources/Analysis_Results/data_summary/write_individual_srat_object_paths/20200717.v1/Seurat_Object_Paths.20200717.v1.tsv")
## input cell type markers
gene2celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/combine_pt_with_emt_markers/20200911.v1/Kidney_Specific_EMT_Genes.20200911.v1.tsv")
## input id meta data table
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")

# specify thresholds ------------------------------------------------------

for (aliquot2process in "CPT0001260013") {
  aliquot_show <- idmetadata_df$Aliquot.snRNA.WU[idmetadata_df$Aliquot.snRNA == aliquot2process]
  
  # input seurat object, and subset to cluster -----------------------------------------------------
  path_srat <- paths_srat_df$Path_box_seurat_object[paths_srat_df$Aliquot == aliquot2process]
  srat <- readRDS(file = path_srat)
  
  # map cell type ------------------------------------------------------------------
  ## subset barcode2celltype
  barcode2celltype_filtered_df <- barcode2celltype_df %>%
    filter(orig.ident == aliquot2process)
  ## identify major cell type in each cluster
  barcode2celltype_filtered_df$id_seurat_cluster <- mapvalues(x = barcode2celltype_filtered_df$individual_barcode, from = rownames(srat@meta.data), to = as.vector(srat@meta.data$seurat_clusters))
  table(barcode2celltype_filtered_df$id_seurat_cluster)
  table(srat@meta.data$seurat_clusters)
  count_bycelltype_bycluster_df <- barcode2celltype_filtered_df %>%
    select(Cell_group3, Cell_group15, id_seurat_cluster) %>%
    table() %>%
    data.frame() %>%
    filter(Freq > 0)
  num_transitionalcells <- count_bycelltype_bycluster_df$Freq[count_bycelltype_bycluster_df$Cell_group15 == "Transitional cells"]
  num_transitionalcells
  num_nephronepithelium <- sum(count_bycelltype_bycluster_df$Freq[count_bycelltype_bycluster_df$Cell_group3 == "Nephron_Epithelium"])
  num_nephronepithelium
}

## examine normal epithelial cells
normalcells_df <- barcode2celltype_df %>%
  filter(Cell_type.shorter == "Normal epithelial cells")

normallikecells_df <- barcode2celltype_df %>%
  filter(Cell_type.shorter == "Normal-like epithelial cells")
