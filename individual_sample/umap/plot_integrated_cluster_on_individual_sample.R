# Yige Wu @WashU March 2020
## for plotting the cluster number from the integrated data to individual sample clusters

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
## input path to individual srat objects
srat_paths <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/individual_sample/write_path_to_seurat_objects_on_box/20200219.v1/Seurat_Object_Paths.20200219.v1.tsv", data.table = F)
## input barcode to integrated cluster number info
barcode2integrated_cluster_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/fetch_data/20200212.v3/30_aliquot_integration.20200212.v3.umap_data.tsv", data.table = F)

# process the barcode in integrated data to the barcode in individual sample--------------------
barcode2integrated_cluster_df <- barcode2integrated_cluster_df %>%
  rename(integrated_barcode = barcode) %>%
  mutate(individual_barcode = str_split_fixed(string = integrated_barcode, pattern = "_", n = 2)[,1])

# plot cell types to umap by each aliquot ---------------------------------
for (snRNA_aliquot_id_tmp in unique(srat_paths$Aliquot)) {
  ## input seurat object
  seurat_obj_path <- srat_paths$Path_seurat_object[srat_paths$Aliquot == snRNA_aliquot_id_tmp]
  seurat_obj <- readRDS(file = seurat_obj_path)
  
  ## get the coordinates for each cluster label
  p <- DimPlot(seurat_obj, reduction = "umap", label = T, label.size	= 5, repel = T)
  label_data_df <- p$layers[[2]]$data
  
  ## get the coordinates for each cell
  plot_data_df <- FetchData(seurat_obj, vars = c("orig.ident", "ident", "UMAP_1", "UMAP_2"))
  plot_data_df$individual_barcode <- rownames(plot_data_df)
  
  ## map barcode to cell type for this aliquot
  ### get the cell types for all the barcodes in this aliquot
  aliquot_barcode2integrated_cluster_df <- barcode2integrated_cluster_df %>%
    filter(orig.ident == snRNA_aliquot_id_tmp)
  ### map barcode to cell type
  plot_data_df$integrated_cluster_num <- mapvalues(x = plot_data_df$individual_barcode, from = aliquot_barcode2integrated_cluster_df$individual_barcode, to = as.vector(aliquot_barcode2integrated_cluster_df$ident))
  
  ## get the case id for this aliquot to show in the title
  case_id_tmp <- srat_paths$Case[srat_paths$Aliquot == snRNA_aliquot_id_tmp]
  
  ## ggplot
  p <- ggplot()
  
  p <- p + geom_point(data = plot_data_df[!is.na(plot_data_df$cell_Type2Show),], 
                      mapping = aes(x = UMAP_1, y = UMAP_2, color=cell_Type2Show), alpha = 1, size = 0.8)
  p <- p + scale_color_manual(values = colors_immune_cell_types, na.value = "grey80")
  p <- p + ggtitle(paste0("Case: ", case_id_tmp, "   Aliquot: ",  snRNA_aliquot_id_tmp), 
                   subtitle = paste0("Mapping of Immune Cell Types"))
  p <- p + geom_text_repel(data = label_data_df, mapping = aes(UMAP_1, UMAP_2, label = ident), size = 5, colour = "white")
  p <- p +
    theme_bw() +
    theme(panel.border = element_blank(), 
          panel.background = element_rect(fill = "black", colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  p <- p + theme(legend.position = "right",
                 legend.text = element_text(size = 9))
  p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                          axis.text.y=element_blank(),axis.ticks=element_blank(),
                          axis.title.x=element_blank(),
                          axis.title.y=element_blank())
  p
  file2write <- paste(dir_out, snRNA_aliquot_id_tmp, ".Immune_CellTypes_Mapped_From_Integrated_Data.", run_id, ".png", sep="")
  png(file2write, width = 1000, height = 800, res = 150)
  print(p)
  dev.off()
  
  p <- p + geom_point(data = plot_data_df[is.na(plot_data_df$cell_Type2Show),],
                      mapping = aes(x = UMAP_1, y = UMAP_2, color=cell_Type2Show), alpha = 0.5, size = 0.3)
  file2write <- paste(dir_out, snRNA_aliquot_id_tmp, ".Immune_CellTypes_Mapped_From_Integrated_Data_and_notmapped.", run_id, ".png", sep="")
  png(file2write, width = 1000, height = 800, res = 150)
  print(p)
  dev.off()
  
  stop("")
}
