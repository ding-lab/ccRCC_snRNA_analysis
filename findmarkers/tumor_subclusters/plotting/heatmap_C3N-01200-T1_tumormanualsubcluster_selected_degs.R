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
## input srat object
srat_paths <- fread(input = "./Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/recluster_nephron_epithelium_cells_in_individual_samples/20200225.v1/Seurat_Object_Paths.Malignant_Nephron_Epithelium20200225.v1.tsv", data.table = F)
## set aliquot to plot
aliquot_tmp <- "CPT0075130004"
## set genes 2 plot
genes2plot <- c("NFIL3","IMPA2","SP3","ST3GAL1","PFKFB4","FBXO31")


# input srat object and edit meta data ------------------------------------
## get barcode to cnv info for this aliquot
barcode2manualsubcluster_aliquot_df <- barcode2manualsubcluster_df %>%
  filter(orig.ident == aliquot_tmp) %>%
  mutate(Name_Cluster = paste0("C", (Id_TumorManualCluster+1)))
unique(barcode2manualsubcluster_aliquot_df$Name_Cluster)
## input seurat object
seurat_obj_path <- srat_paths$Path_seurat_object[srat_paths$Aliquot == aliquot_tmp]
seurat_obj_path
srat <- readRDS(file = seurat_obj_path)
## change meta data
metadata_df <- srat@meta.data
metadata_df$Name_Cluster <- mapvalues(x = rownames(metadata_df), from = barcode2manualsubcluster_aliquot_df$individual_barcode, to = as.vector(barcode2manualsubcluster_aliquot_df$Name_Cluster))
srat@meta.data <- metadata_df
### set the identities to cluster in the meta data
Idents(object = srat) <- "Name_Cluster"


# plot --------------------------------------------------------------------
FeaturePlot(object = srat, features = genes2plot, order = T)
VlnPlot(object = srat, features = genes2plot, group.by = "Name_Cluster")

p <- DoHeatmap(object = srat, features = genes2plot, group.by = "ident")


png(filename = paste0(dir_out, aliquot_tmp, ".", pathway_tmp, ".", run_id, ".png"),width = 1200, height = 300, res = 150)
print(p)
dev.off()





