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
srat_paths$Path_relative <- gsub(x = srat_paths$Path_seurat_object, pattern = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/", replacement = "./")
srat_paths$Path_relative
## input the barcode-manualsubcluster info
barcode2subclusterid_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_corrected_by_individual_sample_inspection/20200811.v1/31Aliquot.Barcode2CellType.20200811.v1.tsv", data.table = F)
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200505.v1/meta_data.20200505.v1.tsv")

# set parameters for findmarkers ------------------------------------------
## set min.pct
logfc.threshold.wilcox <- 0.1
min.pct.wilcox <- 0.1

# change metadata ---------------------------------------------------------
## format the barcode table
barcode2subclusterid_df$id_aliquot_wu <- mapvalues(x = barcode2subclusterid_df$orig.ident, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
barcode2subclusterid_df <- barcode2subclusterid_df %>%
  mutate(id_aliquot_cluster = paste0(id_aliquot_wu, "_C", (Id_TumorManualCluster + 1)))

# for each aliquot, input seurat object and fetch data and write data --------------------
aliquot_wu_tmp <- "C3L-00010-T1"
aliquot_tmp <- idmetadata_df$Aliquot.snRNA[idmetadata_df$Aliquot.snRNA.WU == aliquot_wu_tmp]

## input individual seurat object
seurat_obj_path <- srat_paths$Path_relative[srat_paths$Aliquot == aliquot_tmp]
seurat_obj_path
srat <- readRDS(file = seurat_obj_path)

### change meta data
barcode2subclusterid_aliquot_df <- barcode2subclusterid_df %>%
  filter(id_aliquot_wu == aliquot_wu_tmp) %>%
  filter(individual_barcode %in% rownames(srat@meta.data))
srat@meta.data$id_aliquot_cluster <- mapvalues(x = rownames(srat@meta.data), from = barcode2subclusterid_aliquot_df$individual_barcode, to = as.vector(barcode2subclusterid_aliquot_df$id_aliquot_cluster))

### change ident
Idents(srat) <- "id_aliquot_cluster"

## specify the groups to compare
id_aliquot_cluster.1 <- c("C3L-00010-T1_C1", "C3L-00010-T1_C2")
id_aliquot_cluster.2 <- c("C3L-00010-T1_C3", "C3L-00010-T1_C4")

## find all markers using Wilcox testing
markers_wilcox_df <- FindMarkers(object = srat, test.use = "wilcox", only.pos = F, min.pct = min.pct.wilcox, logfc.threshold = logfc.threshold.wilcox,
                                 ident.1 = id_aliquot_cluster.1, ident.2 = id_aliquot_cluster.2, verbose = T)
markers_wilcox_df$deg_gene_symbol <- rownames(markers_wilcox_df)
markers_wilcox_df$id_aliquot_wu <- aliquot_wu_tmp
markers_wilcox_df$ident.1 <- paste0(id_aliquot_cluster.1, collapse = ";")
markers_wilcox_df$ident.2 <- paste0(id_aliquot_cluster.2, collapse = ";")

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "C3N-00010.C12_vs_C34.Tumormanualsubcluster.FindMarkers.Wilcox.Minpct", min.pct.wilcox, ".Logfc", logfc.threshold.wilcox,".tsv")
write.table(markers_wilcox_df, file = file2write, 
            quote = F, sep = "\t", row.names = F)