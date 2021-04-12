# Yige Wu @WashU Feb 2020

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
## input the paths for reclustered seurat objects
srat_paths <- fread(input = "./Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/recluster_nephron_epithelium_cells_in_individual_samples/20200225.v1/Seurat_Object_Paths.Malignant_Nephron_Epithelium20200225.v1.tsv", data.table = F)
## input the barcode-manualsubcluster info
barcode2subclusterid_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_tumorsubclusterid/map_barcode_with_manual_tumorsubcluster_id/20201130.v1/Barcode2TumorSubclusterId.20201130.v1.tsv", data.table = F)
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20200902.v1/scrublet.run20200902_adj_cutoff.united_outputs.tsv", data.table = F)

## set min.pct
# min.pct.wilcox <- 0.25
min.pct.wilcox <- 0.1
logfc.threshold.wilcox <- 0.1
# for each aliquot, input seurat object and fetch data and write data --------------------
markers_wilcox_df <- NULL
for (aliquot_tmp in unique(srat_paths$Aliquot)) {
  ## input seurat object
  seurat_obj_path <- srat_paths$Path_seurat_object[srat_paths$Aliquot == aliquot_tmp]
  seurat_obj_path
  srat <- readRDS(file = seurat_obj_path)
  ## change ident to manual subcluster
  ### get the barcode-subclusterid for this aliquot
  metadata_new_df <- barcode2subclusterid_df %>%
    filter(orig.ident == aliquot_tmp)
  rownames(metadata_new_df) <- metadata_new_df$barcode
  ### change meta data
  srat@meta.data <- metadata_new_df
  ### change ident
  Idents(srat) <- "manual_cluster_id"
  
  ## find all markers using Wilcox testing
  markers_wilcox_tmp <- FindAllMarkers(object = srat, test.use = "wilcox", only.pos = F, min.pct = min.pct.wilcox, logfc.threshold = logfc.threshold.wilcox, return.thresh = 0.05)
  
  ## add the current DEGs into the super table
  if (nrow(markers_wilcox_tmp) > 0) {
    markers_wilcox_tmp$id_aliquot <- aliquot_tmp
    markers_wilcox_df <- rbind(markers_wilcox_df, markers_wilcox_tmp)
  }
}

# write output ------------------------------------------------------------
# write.table(markers_wilcox_df, file = paste0(dir_out, "Tumormanualsubcluster.FindAllMarkers.Wilcox.Pos.Minpct", min.pct.wilcox, ".Logfc", logfc.threshold.wilcox,".tsv"),
write.table(markers_wilcox_df, file = paste0(dir_out, "Tumormanualsubcluster.FindAllMarkers.Wilcox.Minpct", min.pct.wilcox, ".Logfc", logfc.threshold.wilcox,".tsv"), 
            quote = F, sep = "\t", row.names = F)
