# Yige Wu @WashU Feb 2020
## fetch UMAP coordinates and cluster info for the tumor cell reclustering

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
## input the paths for reclustered seurat objects
srat_paths <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/recluster_nephron_epithelium_cells_in_individual_samples/20200225.v1/Seurat_Object_Paths.Malignant_Nephron_Epithelium20200225.v1.tsv", data.table = F)
## inpu the integrated barocde name
barcode2celltype_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/map_celltype_to_barcode/20200309.v1/30_aliquot_integration.barcode2celltype.20200309.v1.tsv", data.table = F)

# for each aliquot, input seurat object and fetch data and write data --------------------
aliquots_data <- NULL
for (aliquot_tmp in unique(srat_paths$Aliquot)) {
  ## input seurat object
  seurat_obj_path <- srat_paths$Path_seurat_object[srat_paths$Aliquot == aliquot_tmp]
  seurat_obj_path
  seurat_object <- readRDS(file = seurat_obj_path)
  ## fetch data
  aliquot_data <- FetchData(object = seurat_object, vars = c("orig.ident", "ident", "UMAP_1", "UMAP_2"))
  ## add barcode column
  aliquot_data$individual_barcode <- rownames(aliquot_data)
  ## merge data across aliquots
  aliquots_data <- rbind(aliquots_data, aliquot_data)
}

# rename and add integrated barcode ---------------------------------------
aliquots_data <- merge(aliquots_data, barcode2celltype_df %>%
                         select(orig.ident, individual_barcode, integrated_barcode),
                       by = c("orig.ident", "individual_barcode"), all.x = T)
aliquots_data <- aliquots_data %>%
  rename(tumor_exp_subcluster.ident = ident)


# write table -------------------------------------------------------------
write.table(x = aliquots_data, file = paste0(dir_out, "Tumor_Subcluster_UMAP_Data", ".", run_id, ".tsv"), row.names = F, quote = F, sep = "\t")

