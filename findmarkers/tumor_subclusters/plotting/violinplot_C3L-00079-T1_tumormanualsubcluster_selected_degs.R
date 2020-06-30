# Yige Wu @WashU Jun 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input barcode-to-manual grouping 
barcode2manualsubcluster_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_barcode_with_manual_tumorsubcluster_id/20200616.v1/Barcode2TumorSubclusterId.20200616.v1.tsv", data.table = F)
## load meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200427.v1/meta_data.20200427.v1.tsv", data.table = F)
## input paths srat object
srat_paths <- fread(input = "./Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/recluster_nephron_epithelium_cells_in_individual_samples/20200225.v1/Seurat_Object_Paths.Malignant_Nephron_Epithelium20200225.v1.tsv", data.table = F)
## set sample id
id_aliquot_wu_filter <- "C3L-00079-T1"
## input selected deg
genes2plot_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/examine_degs/examine_tumormanualsubcluster_degs_C3L-00079-T1/20200629.v1/C3L-00079-T1.DEG_selected.tsv")

# get ids -----------------------------------------------------------------
## set aliquot to plot
id_aliquot <- idmetadata_df$Aliquot.snRNA[idmetadata_df$Aliquot.snRNA.WU == id_aliquot_wu_filter]
id_aliquot
## set genes 2 plot
genes2plot <- c("HIF1A", unique(genes2plot_df$deg_gene_symbol))

# input srat object and edit meta data ------------------------------------
## get barcode to cnv info for this aliquot
barcode2manualsubcluster_aliquot_df <- barcode2manualsubcluster_df %>%
  filter(orig.ident == id_aliquot) %>%
  mutate(Name_Cluster = paste0("C", (Id_TumorManualCluster+1)))
unique(barcode2manualsubcluster_aliquot_df$Name_Cluster)
## input seurat object
seurat_obj_path <- srat_paths$Path_seurat_object[srat_paths$Aliquot == id_aliquot]
seurat_obj_path
srat <- readRDS(file = seurat_obj_path)
## change meta data
metadata_df <- srat@meta.data
metadata_df$Name_Cluster <- mapvalues(x = rownames(metadata_df), from = barcode2manualsubcluster_aliquot_df$individual_barcode, to = as.vector(barcode2manualsubcluster_aliquot_df$Name_Cluster))
srat@meta.data <- metadata_df
### set the identities to cluster in the meta data
Idents(object = srat) <- "Name_Cluster"

# plot --------------------------------------------------------------------
for (gene_tmp in genes2plot) {
  p <- VlnPlot(object = srat, features = gene_tmp, group.by = "Name_Cluster", idents = c("C1", "C2", "C3", "C4"), pt.size = 0.25)
  png(filename = paste0(dir_out, id_aliquot_wu_filter, ".", gene_tmp, ".png"),width = 600, height = 500, res = 150)
  print(p)
  dev.off()
}






