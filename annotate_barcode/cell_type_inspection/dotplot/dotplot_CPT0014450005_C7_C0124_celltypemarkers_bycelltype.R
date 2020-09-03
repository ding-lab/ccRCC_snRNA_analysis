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
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_to_all_cells/20200720.v1/30AliquotIntegration.Barcode2CellType.TumorManualCluster.20200720.v1.tsv", data.table = F)
## input seurat paths
paths_srat_df <- fread(data.table = F, input = "./Resources/Analysis_Results/data_summary/write_individual_srat_object_paths/20200717.v1/Seurat_Object_Paths.20200717.v1.tsv")
## input cell type markers
gene2celltype_df <- fread(data.table = F, input = "./Resources/Knowledge/Kidney_Markers/Gene2CellType_Tab.20200406.v1.tsv")
## specify the aliquot id
aliquot2process <- "CPT0014450005"
## specify the cluster id
seuratcluster2process <- c(7, 0, 1, 2, 4)

# input seurat object, and subset to cluster -----------------------------------------------------
path_srat <- paths_srat_df$Path_box_seurat_object[paths_srat_df$Aliquot == aliquot2process]
srat <- readRDS(file = path_srat)
srat <- subset(srat, idents = seuratcluster2process)

# map cell type ------------------------------------------------------------------
## subset barcode2celltype
barcode2celltype_filtered_df <- barcode2celltype_df %>%
  filter(orig.ident == aliquot2process)
## map cell type
srat@meta.data$Cell_type.detailed <- mapvalues(x = rownames(srat@meta.data), from = barcode2celltype_filtered_df$individual_barcode, to = as.vector(barcode2celltype_filtered_df$Cell_type.detailed))
Idents(srat) <- "Cell_type.detailed"
table(srat@meta.data$Cell_type.detailed)
srat <- subset(srat, idents = c("Fibroblasts", "Tumor cells"))
srat@meta.data$Id_Cluster_Celltype <- paste0(srat@meta.data$seurat_clusters, "_", srat@meta.data$Cell_type.detailed)
Idents(srat) <- "Id_Cluster_Celltype"

# prepare data ------------------------------------------------------------
## get the genes within the cell type marker table
genes2plot <-  intersect(gene2celltype_df$Gene, srat@assays$RNA@counts@Dimnames[[1]])
genes2plot <- unique(genes2plot)
## get the pct expressed for each gene in each cluster
p <- DotPlot(object = srat, features = genes2plot, col.min = 0)
plot_data <- p$data
## transform the dataframe to matrix to better filter out genes with too low expressin
plot_matrix <- dcast(data = plot_data, formula = features.plot ~ id, value.var = "pct.exp")
plot_matrix %>% head()
## filter for genes that are expressed in >25% of one cluster at least
## replot with the filtered genes plus malignant cell marker genes
malignant_markers <- as.vector(gene2celltype_df$Gene[gene2celltype_df$Cell_Type_Group == "Malignant_Nephron_Epithelium"])
genes2plot_filtered <- as.vector(plot_matrix[rowSums(plot_matrix[,unique(as.vector(plot_data$id))] > 20) >= 1, "features.plot"])
genes2plot_filtered <- c(genes2plot_filtered, 
                         as.vector(plot_matrix[(plot_matrix$features.plot %in% malignant_markers) & (rowSums(plot_matrix[,unique(as.vector(plot_data$id))] > 10) >= 1), "features.plot"]))
genes2plot_filtered <- unique(genes2plot_filtered)


# plot scaled -------------------------------------------------------------
p <- DotPlot(object = srat, features = genes2plot_filtered, col.min = 0)
p$data$gene_cell_type_group <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type_Group)
p$data$gene_cell_type1 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type1)
p$data$gene_cell_type2 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type2)
p$data$gene_cell_type3 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type3)
p$data$gene_cell_type4 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type4)
p$data$cell_type <- str_split_fixed(string = p$data$id, pattern = "_", n = 2)[,2]
p$data$cluster <- str_split_fixed(string = p$data$id, pattern = "_", n = 2)[,1]
p$data$cell_group <- plyr::mapvalues(p$data$cell_type, from = barcode2celltype_filtered_df$Cell_type.detailed, to = barcode2celltype_filtered_df$Cell_group)
cat("###########################################\n")
cat("Dotplot now\n")
p <- p  + RotatedAxis()
p <- p + facet_grid(cell_group~gene_cell_type_group + gene_cell_type1 + gene_cell_type2 + gene_cell_type3 + gene_cell_type4, scales = "free", space = "free", drop = T)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               panel.grid.major = element_line(colour = "grey80"),
               strip.text.x = element_text(angle = 0, vjust = 0.5),
               strip.text.y = element_text(angle = 0, vjust = 0.5),
               axis.text.x = element_text(size = 15, face = "bold"),
               strip.placement = "outside")
file2write <- paste0(dir_out, aliquot2process, ".C", seuratcluster2process, ".CellTypeMarkerExp.Scaled.png")
png(file = file2write, width = 3500, height = 700, res = 150)
print(p)
dev.off()

# plot not scaled -------------------------------------------------------------
p <- DotPlot(object = srat, features = genes2plot_filtered, col.max = 7, scale = F)
p$data$gene_cell_type_group <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type_Group)
p$data$gene_cell_type1 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type1)
p$data$gene_cell_type2 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type2)
p$data$gene_cell_type3 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type3)
p$data$gene_cell_type4 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type4)
p$data$cell_type <- str_split_fixed(string = p$data$id, pattern = "_", n = 2)[,2]
p$data$cluster <- str_split_fixed(string = p$data$id, pattern = "_", n = 2)[,1]
p$data$cell_group <- plyr::mapvalues(p$data$cell_type, from = barcode2celltype_filtered_df$Cell_type.detailed, to = barcode2celltype_filtered_df$Cell_group)
cat("###########################################\n")
cat("Dotplot now\n")
p <- p +scale_color_gradient2(midpoint=median(p$data$avg.exp, na.rm = T), low="blue", mid="white",
                              high="red", space ="Lab" )
p <- p  + RotatedAxis()
p <- p + facet_grid(cell_group~gene_cell_type_group + gene_cell_type1 + gene_cell_type2 + gene_cell_type3 + gene_cell_type4, scales = "free", space = "free", drop = T)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               panel.grid.major = element_line(colour = "grey80"),
               strip.text.x = element_text(angle = 0, vjust = 0.5),
               strip.text.y = element_text(angle = 0, vjust = 0.5),
               axis.text.x = element_text(size = 15, face = "bold"),
               strip.placement = "outside")
file2write <- paste0(dir_out, aliquot2process, ".C", seuratcluster2process, ".CellTypeMarkerExp.NotScaled.png")
png(file = file2write, width = 3500, height = 700, res = 150)
print(p)
dev.off()