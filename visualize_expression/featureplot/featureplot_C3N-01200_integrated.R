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
## input id meta data table
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")
## input cell type markers
gene2celltype_df <- data.frame(Gene = c("MRC1", "CD163", "MSR1",
                                        "CD86", "FCGR3A", "TLR2", "SRSF10", "SIGLEC1"),
                               Gene_Group = c(rep(x = "M2", 3), rep(x = "M1", 5)))
aliquot2process <- idmetadata_df$Aliquot.snRNA[idmetadata_df$Case == "C3N-01200"]

# input seurat object, and subset to cluster -----------------------------------------------------
## input srat object
srat <- readRDS(file = "./Resources/snRNA_Processed_Data/Merged_Seurat_Objects/C3N-01200.Tumor_Segments.Merged.20200319.v1.RDS")

# annotate cell tyep ------------------------------------------------------
barcode2celltype_aliquot_df <- barcode2celltype_df %>%
  filter(orig.ident == aliquot2process)
srat@meta.data$Cell_type.shorter <- mapvalues(x = rownames(srat@meta.data), from = barcode2celltype_aliquot_df$individual_barcode, to = as.vector(barcode2celltype_aliquot_df$Cell_type.shorter))
srat@meta.data$Cell_group.detailed <- mapvalues(x = rownames(srat@meta.data), from = barcode2celltype_aliquot_df$individual_barcode, to = as.vector(barcode2celltype_aliquot_df$Cell_group.detailed))
srat@meta.data$easy_id <- mapvalues(x = srat@meta.data$orig.ident, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))

# specify genes to plot ---------------------------------------------------
genes_plot <- gene2celltype_df$Gene

# plot by gene ------------------------------------------------------------
for (gene_plot in genes_plot) {
  p <- FeaturePlot(object = srat, features = gene_plot, order = T, min.cutoff = "q1", max.cutoff = "q99", cols = c("yellow", "red"), split.by = "easy_id")
  ## write output
  file2write <- paste0(dir_out, gene_plot, ".png")
  png(file2write, width = 2400, height = 800, res = 150)
  print(p)
  dev.off()
}

