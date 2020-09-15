# Yige Wu @WashU Apr 2020

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
## input cell type per barcode table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups/20200910.v2/31Aliquot.Barcode2CellType.20200910.v2.tsv", data.table = F)
## input cell type markers
gene2celltype_df <- fread(data.table = F, input = "./Resources/Knowledge/Kidney_Markers/Gene2CellType_Tab.20200831.v1.tsv")
## input id meta data table
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")
aliquot_show <- "C3N-01200 Merged"

# input seurat object, add cell type -----------------------------------------------------
srat <- readRDS(file = "./Resources/snRNA_Processed_Data/Merged_Seurat_Objects/C3N-01200.Tumor_Segments.Merged.20200319.v1.RDS")
srat@meta.data$barcode_individual <- str_split_fixed(string = rownames(srat@meta.data), pattern = "_", n = 3)[,1]
srat@meta.data$id_cell <- paste0(srat@meta.data$orig.ident, "_", srat@meta.data$barcode_individual)
## add id_cell to cell type table
barcode2celltype_df <- barcode2celltype_df %>%
  filter(orig.ident %in% unique(srat@meta.data$orig.ident)) %>%
  mutate(id_cell = paste0(orig.ident, '_', individual_barcode))
## map cell type
srat@meta.data$Cell_group3 <- mapvalues(x = srat@meta.data$id_cell, from = barcode2celltype_df$id_cell, to = barcode2celltype_df$Cell_group3)
srat@meta.data$Cell_type.shorter <- mapvalues(x = srat@meta.data$id_cell, from = barcode2celltype_df$id_cell, to = barcode2celltype_df$Cell_type.shorter)
## add 
srat@meta.data$Aliquot_WU <- mapvalues(x = srat@meta.data$orig.ident, from = idmetadata_df$Aliquot.snRNA, to = idmetadata_df$Aliquot.snRNA.WU)
srat@meta.data$Id_PlotGroup <- paste0("C", srat@meta.data$seurat_clusters, "_", srat@meta.data$Cell_type.shorter, '_', srat@meta.data$Aliquot_WU)
Idents(srat) <- "Cell_group3"

for (cellgroup_plot in c("Nephron_Epithelium")) {
  ## subset
  srat_plot <- subset(x = srat, idents = cellgroup_plot)
  dim(srat_plot)
  unique(srat_plot@meta.data$Cell_type.shorter)
  ## count plot group
  countcells_plotgroup <- data.frame(table(srat_plot@meta.data$Id_PlotGroup))
  ## subset again
  countcells_plotgroup_keep <- countcells_plotgroup$Var1[countcells_plotgroup$Freq >= 10]
  Idents(srat_plot) <- "Id_PlotGroup"
  srat_plot <- subset(x = srat_plot, idents = countcells_plotgroup_keep)
  dim(srat_plot)

  # prepare data ------------------------------------------------------------
  ## get the genes within the cell type marker table
  genes2plot <-  intersect(gene2celltype_df$Gene, srat_plot@assays$RNA@counts@Dimnames[[1]])
  genes2plot <- unique(genes2plot)
  ## get the pct expressed for each gene in each cluster
  p <- DotPlot(object = srat_plot, features = genes2plot, col.min = 0, assay = "RNA")
  plot_data <- p$data
  ## transform the dataframe to matrix to better filter out genes with too low expressin
  plot_matrix <- dcast(data = plot_data, formula = features.plot ~ id, value.var = "pct.exp")
  plot_matrix %>% head()
  ## filter for genes that are expressed in >25% of one cluster at least
  ## replot with the filtered genes plus malignant cell marker genes
  malignant_markers <- as.vector(gene2celltype_df$Gene[gene2celltype_df$Cell_Type_Group == "Malignant_Nephron_Epithelium"])
  genes2plot_filtered <- as.vector(plot_matrix[rowSums(plot_matrix[,unique(as.vector(plot_data$id))] > 50) >= 1, "features.plot"])
  genes2plot_filtered <- c(genes2plot_filtered,
                           as.vector(plot_matrix[(plot_matrix$features.plot %in% malignant_markers) & (rowSums(plot_matrix[,unique(as.vector(plot_data$id))] > 50) >= 1), "features.plot"]))
  genes2plot_filtered <- unique(genes2plot_filtered)

  # plot scaled -------------------------------------------------------------
  p <- DotPlot(object = srat_plot, features = genes2plot_filtered, col.min = 0, assay = "RNA")
  p$data$gene_cell_type_group <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type_Group)
  p$data$gene_cell_type1 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type1)
  p$data$gene_cell_type2 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type2)
  p$data$gene_cell_type3 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type3)
  p$data$gene_cell_type4 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type4)
  p$data$cell_type <- str_split_fixed(string = p$data$id, pattern = "_", n = 2)[,2]
  p$data$cluster <- str_split_fixed(string = p$data$id, pattern = "_", n = 2)[,1]
  p$data$cell_group <- plyr::mapvalues(p$data$cell_type, from = barcode2celltype_filtered_df$Cell_type.detailed, to = barcode2celltype_filtered_df$Cell_group.detailed)
  cat("###########################################\n")
  cat("Dotplot now\n")
  p <- p  + RotatedAxis()
  p <- p + facet_grid(cluster~gene_cell_type_group + gene_cell_type1 + gene_cell_type2 + gene_cell_type3 + gene_cell_type4, scales = "free", space = "free", drop = T)
  p <- p + theme(panel.spacing = unit(0, "lines"),
                 strip.background = element_blank(),
                 panel.border = element_rect(colour = "black"),
                 panel.grid.major = element_line(colour = "grey80"),
                 strip.text.x = element_text(angle = 0, vjust = 0.5),
                 strip.text.y = element_text(angle = 0, vjust = 0.5),
                 axis.text.x = element_text(size = 15, face = "bold"),
                 strip.placement = "outside")
  p <- p + ggtitle(paste0(aliquot_show, " Expression of Cell Type Marker Genes"))
  file2write <- paste0(dir_out, aliquot_show, ".CellTypeMarkerExp.Scaled.png")
  png(file = file2write, width = 5000, height = 1500, res = 150)
  print(p)
  dev.off()

  # plot not scaled -------------------------------------------------------------
  p <- DotPlot(object = srat_plot, features = genes2plot_filtered, col.max = 7, scale = F, assay = "RNA")
  p$data$gene_cell_type_group <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type_Group)
  p$data$gene_cell_type1 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type1)
  p$data$gene_cell_type2 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type2)
  p$data$gene_cell_type3 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type3)
  p$data$gene_cell_type4 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type4)
  p$data$cell_type <- str_split_fixed(string = p$data$id, pattern = "_", n = 2)[,2]
  p$data$cluster <- str_split_fixed(string = p$data$id, pattern = "_", n = 2)[,1]
  p$data$cell_group <- plyr::mapvalues(p$data$cell_type, from = barcode2celltype_filtered_df$Cell_type.detailed, to = barcode2celltype_filtered_df$Cell_group.detailed)
  cat("###########################################\n")
  cat("Dotplot now\n")
  p <- p +scale_color_gradient2(midpoint=median(p$data$avg.exp, na.rm = T), low="blue", mid="white",
                                high="red", space ="Lab" )
  p <- p  + RotatedAxis()
  p <- p + facet_grid(cluster~gene_cell_type_group + gene_cell_type1 + gene_cell_type2 + gene_cell_type3 + gene_cell_type4, scales = "free", space = "free", drop = T)
  p <- p + theme(panel.spacing = unit(0, "lines"),
                 strip.background = element_blank(),
                 panel.border = element_rect(colour = "black"),
                 panel.grid.major = element_line(colour = "grey80"),
                 strip.text.x = element_text(angle = 0, vjust = 0.5),
                 strip.text.y = element_text(angle = 0, vjust = 0.5),
                 axis.text.x = element_text(size = 15, face = "bold"),
                 strip.placement = "outside")
  p <- p + ggtitle(paste0(aliquot_show, " Expression of Cell Type Marker Genes"))
  file2write <- paste0(dir_out, aliquot_show, ".CellTypeMarkerExp.NotScaled.png")
  png(file = file2write, width = 5000, height = 1500, res = 150)
  print(p)
  dev.off()
  
}

