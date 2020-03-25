# Yige Wu @WashU March 2020
## for each individual sample tumor cell reclustered, plot UMAP

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/aes.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input seurat object paths
srat_paths <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/recluster_nephron_epithelium_cells_in_individual_samples/20200225.v1/Seurat_Object_Paths.Malignant_Nephron_Epithelium20200225.v1.tsv", data.table = F)

# plot by each aliquot ----------------------------------------------------
for (aliquot_tmp in srat_paths$Aliquot) {
  ## input srat object
  srat_path <- srat_paths$Path_seurat_object[srat_paths$Aliquot == aliquot_tmp]
  srat <- readRDS(file = srat_path)
  
  ## make color for each cluster
  uniq_cluster_colors <- Polychrome::dark.colors(n = length(unique(srat@meta.data$seurat_clusters)))
  names(uniq_cluster_colors) <- unique(srat@meta.data$seurat_clusters)
  
  ## make plot
  p <- DimPlot(object = srat)
  p <- p + ggtitle(label = paste0(aliquot_tmp, " Tumor-Cell-Only Clustering"))
  p <- p + scale_color_manual(values = uniq_cluster_colors)
  p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                          axis.text.y=element_blank(),axis.ticks=element_blank(),
                          axis.title.x=element_blank(),
                          axis.title.y=element_blank())
  ## save plot
  file2write <- paste0(dir_out, aliquot_tmp, ".TumorCellOnlyClustering.", run_id, ".png")
  png(filename = file2write, width = 1000, height = 900, res = 150)
  print(p)
  dev.off()
}