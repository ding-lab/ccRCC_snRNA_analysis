# Yige Wu @WashU Oct 2019
## for generating the barcode annotation files (assigned to one random group) for 10Xmapping
## only use barcodes of after QC barcodes so that 10Xmapping won't run forever

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
## input seurat object master list
srat_paths_df <- fread(input = "./Resources/Analysis_Results/data_summary/write_individual_srat_object_paths/20210729.v1/Seurat_Object_Paths.20210729.v1.tsv", data.table = F)

# input seurat object -----------------------------------------------------
for (i in 1:7) {
# for (i in 1:nrow(srat_paths_df)) {
  snRNA_aliquot_id_tmp <- srat_paths_df$Aliquot[i]
  # facs_tmp <- srat_paths_df$FACS[i]
  path_seurat_obj_tmp <- srat_paths_df$Path_box_seurat_object[i]
  seurat_obj_tmp <- readRDS(file = path_seurat_obj_tmp)
  anno_tab_tmp <- seurat_obj_tmp@meta.data
  anno_tab_tmp$barcode <- rownames(anno_tab_tmp)
  anno_tab_tmp <- anno_tab_tmp %>%
    select(barcode) %>%
    mutate(random_group = "0")
  # write.table(x = anno_tab_tmp, file = paste0(dir_out, snRNA_aliquot_id_tmp, facs_tmp, "_AfterQC_Barcodes.tsv"), quote = F, row.names = F, sep = "\t", col.names = F)
  write.table(x = anno_tab_tmp, file = paste0(dir_out, snRNA_aliquot_id_tmp, "_AfterQC_Barcodes.tsv"), quote = F, row.names = F, sep = "\t", col.names = F)
}


