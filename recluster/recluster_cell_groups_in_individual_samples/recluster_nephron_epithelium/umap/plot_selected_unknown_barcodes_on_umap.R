# Yige Wu @WashU March 2020
## plot the unknown barcode from Song to the umap

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
srat_paths <- srat_paths %>%
  filter(Aliquot %in% c("CPT0075140002", "CPT0075120002"))
## input the unknown barcodes
unknown_barcodes_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/Multi-Segment_Patients/C3N-01200.Tumor_Segments.Merged.20200211.v1.cluster7.tsv", data.table = F)
### preprocess
unknown_barcodes_df <- unknown_barcodes_df %>%
  mutate(individual_barcode = str_split_fixed(string = bc, pattern = "_", n = 2)[,1])
# plot by each aliquot ----------------------------------------------------
for (aliquot_tmp in srat_paths$Aliquot) {
  ## input srat object
  srat_path <- srat_paths$Path_seurat_object[srat_paths$Aliquot == aliquot_tmp]
  srat <- readRDS(file = srat_path)
  
  ## make color for each cluster
  uniq_cluster_colors <- Polychrome::dark.colors(n = length(unique(srat@meta.data$seurat_clusters)))
  names(uniq_cluster_colors) <- unique(srat@meta.data$seurat_clusters)
  
  ## make plot data
  p <- DimPlot(object = srat)
  plot_data_df <- p$data
  plot_data_df$barcode <- rownames(plot_data_df)
  plot_data_df <- plot_data_df %>%
    mutate(is_unknown = (barcode %in% unknown_barcodes_df$individual_barcode[unknown_barcodes_df$Sample == aliquot_tmp]))
  
  ## make plot
  p <- ggplot()
  p <- p + geom_point(data = plot_data_df, mapping = aes(x = UMAP_1, y = UMAP_2, color = is_unknown))
  p <- p + ggtitle(label = paste0(aliquot_tmp, " Tumor-Cell-Only Clustering"))
  p <- p + scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey50"))
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