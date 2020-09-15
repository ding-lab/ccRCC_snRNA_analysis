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
## input seurat paths
paths_srat_df <- fread(data.table = F, input = "./Resources/Analysis_Results/data_summary/write_individual_srat_object_paths/20200717.v1/Seurat_Object_Paths.20200717.v1.tsv")
## input cell type markers
gene2celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/combine_pt_with_emt_markers/20200911.v1/Kidney_Specific_EMT_Genes.20200911.v1.tsv")
## input id meta data table
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")

# specify thresholds ------------------------------------------------------
## filter for genes that are expressed in >25% of one cluster at least
pct_thres <- 50
avgexp_thres <- 1

for (aliquot2process in "CPT0001260013") {
  aliquot_show <- idmetadata_df$Aliquot.snRNA.WU[idmetadata_df$Aliquot.snRNA == aliquot2process]
  
  # input seurat object, and subset to cluster -----------------------------------------------------
  path_srat <- paths_srat_df$Path_box_seurat_object[paths_srat_df$Aliquot == aliquot2process]
  srat <- readRDS(file = path_srat)
  
  # map cell type ------------------------------------------------------------------
  ## subset barcode2celltype
  barcode2celltype_filtered_df <- barcode2celltype_df %>%
    filter(orig.ident == aliquot2process)
  ## identify major cell type in each cluster
  barcode2celltype_filtered_df$id_seurat_cluster <- mapvalues(x = barcode2celltype_filtered_df$individual_barcode, from = rownames(srat@meta.data), to = as.vector(srat@meta.data$seurat_clusters))
  table(barcode2celltype_filtered_df$id_seurat_cluster)
  table(srat@meta.data$seurat_clusters)
  count_bycelltype_bycluster_df <- barcode2celltype_filtered_df %>%
    select(Cell_group3, Cell_group14, id_seurat_cluster) %>%
    table() %>%
    data.frame() %>%
    filter(Freq > 0) %>%
    mutate(Id_Cluster_CellType = paste0("C", id_seurat_cluster, "_", Cell_group14)) %>%
    mutate(Keep = (Freq >= 30 & (Cell_group3 %in% c("Nephron_Epithelium") | Cell_group14 %in% c("Fibroblasts", "Myofibroblasts"))))
  
  barcode2celltype_filtered_df <- barcode2celltype_filtered_df %>%
    mutate(Id_Cluster_CellType = paste0("C", id_seurat_cluster, "_", Cell_group14))
  
  barcode2celltype_filtered_df$Keep <- mapvalues(x = barcode2celltype_filtered_df$Id_Cluster_CellType, from = count_bycelltype_bycluster_df$Id_Cluster_CellType, to = as.vector(count_bycelltype_bycluster_df$Keep))
  ## subset seurat object
  srat <- subset(srat, cells = barcode2celltype_filtered_df$individual_barcode[barcode2celltype_filtered_df$Keep == "TRUE"])
  dim(srat)
  srat@meta.data$Cell_group14 <- mapvalues(x = rownames(srat@meta.data), from = barcode2celltype_filtered_df$individual_barcode, to = barcode2celltype_filtered_df$Cell_group14)
  Idents(srat) <- "Cell_group14"
  
  # prepare data ------------------------------------------------------------
  ## get the genes within the cell type marker table
  genes2plot <-  intersect(gene2celltype_df$Gene, srat@assays$RNA@counts@Dimnames[[1]])
  genes2plot <- unique(genes2plot)
  ## get the pct expressed for each gene in each cluster
  p <- DotPlot(object = srat, features = genes2plot, col.min = 0)
  expdata_df <- p$data
  ## filter genes based on the percentage expressed
  pct_matrix <- dcast(data = expdata_df, formula = features.plot ~ id, value.var = "pct.exp")
  genes_pct_filtered <- as.vector(pct_matrix[rowSums(pct_matrix[,unique(as.vector(expdata_df$id))] > pct_thres) >= 1, "features.plot"])
  ## filter genes based on the average expression
  avgexp_matrix <- dcast(data = expdata_df, formula = features.plot ~ id, value.var = "avg.exp")
  genes_exp_filtered <- as.vector(avgexp_matrix[rowSums(avgexp_matrix[,unique(as.vector(expdata_df$id))] > avgexp_thres) >= 1, "features.plot"])
  ## intersect
  genes2plot_filtered <- intersect(unique(genes_exp_filtered), unique(genes_pct_filtered))
  
  # plot scaled -------------------------------------------------------------
  p <- DotPlot(object = srat, features = genes2plot_filtered, col.min = 0)
  p$data$gene_group1 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Gene_Group1)
  p$data$gene_group2 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Gene_Group2)
  p$data$cell_type <- str_split_fixed(string = p$data$id, pattern = "_", n = 2)[,2]
  p$data$cluster <- str_split_fixed(string = p$data$id, pattern = "_", n = 2)[,1]
  p$data$cell_group <- plyr::mapvalues(p$data$cell_type, from = barcode2celltype_filtered_df$Cell_group14, to = barcode2celltype_filtered_df$Cell_group3)
  cat("###########################################\n")
  cat("Dotplot now\n")
  p <- p  + RotatedAxis()
  p <- p + facet_grid(cell_group~gene_group1 + gene_group2, scales = "free", space = "free", drop = T)
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
  png(file = file2write, width = 3000, height = 600, res = 150)
  print(p)
  dev.off()
  
  # plot not scaled -------------------------------------------------------------
  plotdata_df <- expdata_df %>%
    filter(features.plot %in% genes2plot_filtered)
  expvalue_top <- quantile(x = plotdata_df$avg.exp, probs = 0.95)
  plotdata_df <- plotdata_df %>%
    mutate(expvalue_plot = ifelse(avg.exp >= expvalue_top, expvalue_top, avg.exp))
  
  ## add facet
  plotdata_df$gene_group1 <- plyr::mapvalues(plotdata_df$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Gene_Group1)
  plotdata_df$gene_group2 <- plyr::mapvalues(plotdata_df$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Gene_Group2)
  
  p <- ggplot()
  p <- p + geom_point(data = plotdata_df, mapping = aes(x = features.plot, y = id, color = expvalue_plot, size = pct.exp), shape = 16)
  # p <- p + scale_fill_continuous(type = "viridis")
  p <- p + scale_color_gradientn(colours = rev(RColorBrewer::brewer.pal(n = 9, name = "Spectral")[1:5]))
  p <- p + scale_size_continuous(range = c(0, 8), name="% Expressed")
  # p <- p +scale_color_gradient2(midpoint=median(plotdata_df$avg.exp, na.rm = T), low="blue", mid="white",
                                # high="red", space ="Lab" )
  # p <- p + scale_colour_gradientn(colours = myPalette(100), limits=c(0, 2.5))
  p <- p + facet_grid(.~gene_group1 + gene_group2, scales = "free", space = "free", drop = T)
  p <- p + theme(axis.text.x = element_text(angle = 90, face = "bold", size = 10))
  p <- p + theme(axis.text.y = element_text( face = "bold", size = 12))
  p <- p + theme(panel.spacing = unit(0, "lines"), panel.grid.major = element_line(colour = "grey80"), 
                 panel.border = element_rect(color = "black", fill = NA, size = 0.5),
                 panel.background = element_blank())
  p <- p + theme(strip.background = element_rect(color = "black", fill = NA, size = 0.5), 
                 strip.text.x = element_text(angle = 90, vjust = 0.5), strip.text.y = element_text(angle = 0, vjust = 0.5))
  p <- p + theme(axis.title = element_blank())
  p <- p + ggtitle(paste0(aliquot_show, " Expression of Cell Type Marker Genes"))
  p <- p + labs(colour = "Expression value")
  file2write <- paste0(dir_out, aliquot_show, ".CellTypeMarkerExp.NotScaled.png")
  png(file = file2write, width = 1000, height = 700, res = 150)
  print(p)
  dev.off()
  
  file2write <- paste0(dir_out, aliquot_show, ".CellTypeMarkerExp.NotScaled.legend.png")
  png(file = file2write, width = 1000, height = 1000, res = 150)
  print(p)
  dev.off()
}
