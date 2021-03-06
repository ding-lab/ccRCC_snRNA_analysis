# Yige Wu @WashU Aug 2020

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
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_corrected_by_individual_sample_inspection/20200828.v1/31Aliquot.Barcode2CellType.20200828.v1.tsv", data.table = F)
## input seurat paths
paths_srat_df <- fread(data.table = F, input = "./Resources/Analysis_Results/data_summary/write_individual_srat_object_paths/20200717.v1/Seurat_Object_Paths.20200717.v1.tsv")
## input cell type per barcode table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_corrected_by_individual_sample_inspection/20200828.v1/31Aliquot.Barcode2CellType.20200828.v1.tsv", data.table = F)
## specify the aliquot id
aliquot2process <- "CPT0078510004"

# input seurat object, and subset to cluster -----------------------------------------------------
path_srat <- paths_srat_df$Path_box_seurat_object[paths_srat_df$Aliquot == aliquot2process]
srat <- readRDS(file = path_srat)

# annotate cell tyep ------------------------------------------------------
barcode2celltype_aliquot_df <- barcode2celltype_df %>%
  filter(orig.ident == aliquot2process)
srat@meta.data$Cell_type.shorter <- mapvalues(x = rownames(srat@meta.data), from = barcode2celltype_aliquot_df$individual_barcode, to = as.vector(barcode2celltype_aliquot_df$Cell_type.shorter))
srat@meta.data$Cell_group.detailed <- mapvalues(x = rownames(srat@meta.data), from = barcode2celltype_aliquot_df$individual_barcode, to = as.vector(barcode2celltype_aliquot_df$Cell_group.detailed))

# specify genes to plot ---------------------------------------------------
genes_plot <- c("FN1", "PDGFRB", "LRP2", "PAX8")

# plot by gene ------------------------------------------------------------
for (gene_plot in genes_plot) {
  p <- FeaturePlot(object = srat, features = gene_plot, order = T, min.cutoff = "q1", max.cutoff = "q99", cols = c("yellow", "red"))
  ## write output
  file2write <- paste0(dir_out, gene_plot, ".png")
  png(file2write, width = 1000, height = 1000, res = 150)
  print(p)
  dev.off()
  
  p <- FeaturePlot(object = srat, features = gene_plot, order = T, min.cutoff = "q1", max.cutoff = "q99", cols = c("yellow", "red"), 
                   split.by = "Cell_group.detailed")
  ## write output
  file2write <- paste0(dir_out, gene_plot, ".splitbycellgroup.png")
  png(file2write, width = 5000, height = 1000, res = 150)
  print(p)
  dev.off()
}

