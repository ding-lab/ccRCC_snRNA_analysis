# Yige Wu @WashU Feb 2020
## running on local
## for plotting the dotplot the immune marker genes for the immune reclustered for individual samples

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
## input gene to cell type info
gene2celltype_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Kidney_Markers/Gene2CellType_Tab.20200220.v1.tsv", data.table = F)

# for each aliquot, input immune reclustered object, plot dotplot -----------------------
## get aliquots to input
dir_srat_parent <- "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_immune/recluster_immune_in_individual_samples/20200211.v1/"
aliquots2process <- list.files(path = dir_srat_parent)
## loop
aliquot_tmp <- "CPT0000870003"
for (aliquot_tmp in aliquots2process) {
  ### get the path to the seurat object
  dir_srat <- paste0(dir_srat_parent, aliquot_tmp, "/")
  filename_srat <- list.files(path = dir_srat)
  path_srat <- paste0(dir_srat, filename_srat)
  
  ### input seurat object
  srat <- readRDS(file = path_srat)
  
  ### set genes to plot
  #### get the genes within the cell type marker table
  genes2plot <-  intersect(gene2celltype_df$Gene, srat@assays$RNA@counts@Dimnames[[1]])
  genes2plot <- unique(genes2plot)
  #### get the pct expressed for each gene in each cluster
  p <- DotPlot(object = srat, features = genes2plot, col.min = 0)
  plot_data <- p$data
  #### transform the dataframe to matrix to better filter out genes with too low expressin
  plot_matrix <- dcast(data = plot_data, formula = features.plot ~ id, value.var = "pct.exp")
  plot_matrix %>% head()
  #### filter for genes that are expressed in >25% of one cluster at least
  genes2plot_filtered <- as.vector(plot_matrix[rowSums(plot_matrix[,unique(as.vector(plot_data$id))] > 10) >= 1, "features.plot"])
  
  ### dotplot
  p <- DotPlot(object = srat, features = genes2plot_filtered, col.min = 0)
  p$data$gene_cell_type_group <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type_Group)
  p$data$gene_cell_type1 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type1)
  p$data$gene_cell_type2 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type2)
  p$data$gene_cell_type3 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type3)
  p$data$gene_cell_type4 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type4)
  p <- p  + RotatedAxis()
  p <- p + facet_grid(.~gene_cell_type_group + gene_cell_type1 + gene_cell_type2 + gene_cell_type3 + gene_cell_type4, scales = "free", space = "free", drop = T)
  p <- p + theme(panel.spacing = unit(0, "lines"),
                 strip.background = element_blank(),
                 panel.border = element_rect(colour = "black"),
                 panel.grid.major = element_line(colour = "grey50"),
                 strip.text.x = element_text(angle = 0, vjust = 0.5),
                 axis.text.x = element_text(size = 15, face = "bold"),
                 strip.placement = "outside")
  
  ### save plot
  file2write <- paste0(dir_out, filename_srat, ".dotplot.",run_id, ".png")
  png(file = file2write, width = 4000, height = 1200, res = 150)
  print(p)
  dev.off()
}
