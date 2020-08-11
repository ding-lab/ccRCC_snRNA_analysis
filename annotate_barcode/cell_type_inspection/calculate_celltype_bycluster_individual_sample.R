# Yige Wu @WashU Aug 2020

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
## input cell tpe by barcode
barcode2celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_to_all_cells/20200720.v1/30AliquotIntegration.Barcode2CellType.TumorManualCluster.20200720.v1.tsv")
## input feteched data
barcode2metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/data_summary/fetch_data/fetch_data_by_individual_sample/20200717.v1/Barcode2MetaData.20200717.v1.tsv")

# merge info ------------------------------------------
## merge meta data with cell type
barcode_info_merged_df <- merge(barcode2metadata_df, barcode2celltype_df, by.x = c("aliquot", "individual_barcode"), by.y = c("orig.ident", "individual_barcode"), all.x = T)

# summarize by cluster by cell type shorter by sample ---------------------
summary_bycelltypeshorter_df <- barcode_info_merged_df %>%
  select(aliquot, seurat_cluster_id, Cell_type.shorter) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq >= 50) %>%
  arrange(aliquot, seurat_cluster_id, desc(Freq))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Cells_BySampleByClusterByCellTypeShorter.Over50.", run_id, ".tsv")
write.table(x = summary_bycelltypeshorter_df, file = file2write, quote = F, sep = "\t", row.names = F)

