# Yige Wu @WashU Feb 2020
## for plotting the marker genes for integrated object, showing cell of origin

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
## input seurat processing info to input individual seurat object later
seurat_summary2process <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/individual_cluster/write_path_to_seurat_objects_on_box/20200219.v1/Seurat_Object_Paths.20200219.v1.tsv", data.table = F)
## input cell type assignment for the all immune cells reclustered
cellgroup_integrated_cluster2celltype_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/integration.immune.cluster2celltype.20200214.v1.csv", data.table = F)
## input barcode to cluster info for the all immune cells reclustered
cellgroup_integrated_barcode2cluster_df <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_integrated_data/fetch_data/20200214.v1/integration.202002012.v3.immune_reclustered.20200213.v2.umap_data.20200214.v1.tsv", data.table = F)
## set color for different cell types
colors_immune_cell_types <- c("Macrophages M1" = "#ffff99", ## yellow
                              "Macrophages M1&M2" = "#fdbf6f", ## light orange
                              "Macrophages M2" = "#ff7f00", ## dark orange
                              "cDC1" = "#33a02c", ## dark green
                              "Myeloid lineage immune cells" = "#b2df8a", ## light green
                              "CD4+ T-cells" = "#1f78b4", ## dark blue
                              "CD8+ T-cells" = "#6a3d9a", ## dark purple
                              "NK-cells" = "#cab2d6", ## light purple
                              "CD8+ T-cells & CD4+ T-cells" =  "#a6cee3", ## light blue
                              "B-cells" = "#fb9a99", ## pink
                              "Plasma cells" = "#e31a1c", ## red
                              "Unknown" = "grey50")




# map barcode to cell type ------------------------------------------------
## merge barcode to cluster info with cluster to cell type info
cellgroup_integrated_barcode2celltype_df <- merge(cellgroup_integrated_barcode2cluster_df, cellgroup_integrated_cluster2celltype_df,
                                                  by.x = c("ident"), by.y = c("Cluster"), all.x = T)
## select the cell type to show
celltypes2show <- as.vector(cellgroup_integrated_barcode2celltype_df$Most_Enriched_Cell_Type4)
celltypes2show[celltypes2show == ""] <- as.vector(cellgroup_integrated_barcode2celltype_df$Most_Enriched_Cell_Type3[celltypes2show == ""])
celltypes2show[celltypes2show == ""] <- as.vector(cellgroup_integrated_barcode2celltype_df$Most_Enriched_Cell_Type2[celltypes2show == ""])
celltypes2show[celltypes2show == ""] <- as.vector(cellgroup_integrated_barcode2celltype_df$Most_Enriched_Cell_Type1[celltypes2show == ""])
celltypes2show[celltypes2show == ""] <- as.vector(cellgroup_integrated_barcode2celltype_df$Most_Enriched_Cell_Group[celltypes2show == ""])
cellgroup_integrated_barcode2celltype_df$Cell_Type2Show <- celltypes2show
## format the barcode to merge with the barcode from individual object later
cellgroup_integrated_barcode2celltype_df <- cellgroup_integrated_barcode2celltype_df %>%
  mutate(individual_barcode = str_split_fixed(string = barcode, pattern = "_", n = 2)[,1])

# plot cell types to umap by each aliquot ---------------------------------
for (snRNA_aliquot_id_tmp in unique(seurat_summary2process$Aliquot)) {
  ## input seurat object
  seurat_obj_path <- seurat_summary2process$Path_seurat_object[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
  seurat_obj <- readRDS(file = seurat_obj_path)
  
  ## get the coordinates for each cluster label
  p <- DimPlot(seurat_obj, reduction = "umap", label = T, label.size	= 5, repel = T)
  label_data_df <- p$layers[[2]]$data
  
  ## get the coordinates for each cell
  plot_data_df <- FetchData(seurat_obj, vars = c("orig.ident", "ident", "UMAP_1", "UMAP_2"))
  plot_data_df$barcode <- rownames(plot_data_df)
  
  ## map barcode to cell type for this aliquot
  ### get the cell types for all the barcodes in this aliquot
  aliquot_barcode2celltype_df <- cellgroup_integrated_barcode2celltype_df %>%
    filter(orig.ident == snRNA_aliquot_id_tmp)
  ### map barcode to cell type
  celltype2show <- mapvalues(x = plot_data_df$barcode, from = aliquot_barcode2celltype_df$individual_barcode, to = aliquot_barcode2celltype_df$Cell_Type2Show)
  celltype2show[celltype2show == plot_data_df$barcode] <- NA
  plot_data_df$cell_Type2Show <- celltype2show
  
  ## get the case id for this aliquot to show in the title
  case_id_tmp <- seurat_summary2process$Case[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
  
  
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
}





