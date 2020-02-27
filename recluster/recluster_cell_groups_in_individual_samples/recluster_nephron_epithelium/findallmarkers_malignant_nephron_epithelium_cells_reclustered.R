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

# for each aliquot, input seurat object and fetch data and write data --------------------
for (aliquot_tmp in unique(srat_paths$Aliquot)) {
  ## input seurat object
  seurat_obj_path <- srat_paths$Path_seurat_object[srat_paths$Aliquot == aliquot_tmp]
  seurat_obj_path
  srat <- readRDS(file = seurat_obj_path)
  
  ## find all markers using ROC testing
  markers_roc <- FindAllMarkers(object = srat, test.use = "roc", only.pos = T, return.thresh = 0.5)
  markers_roc <- markers_roc %>%
    filter(power > 0)
  write.table(markers_roc, file = paste0(dir_out, aliquot_tmp, ".FindAllMarkers.ROC.Pos.tsv"), quote = F, sep = "\t", row.names = F)
  
  ## find all markers using Wilcox testing
  markers_wilcox <- FindAllMarkers(object = srat, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
  write.table(markers_wilcox, file = paste0(dir_out, aliquot_tmp, ".FindAllMarkers.Wilcox.Pos.tsv"), quote = F, sep = "\t", row.names = F)
}

