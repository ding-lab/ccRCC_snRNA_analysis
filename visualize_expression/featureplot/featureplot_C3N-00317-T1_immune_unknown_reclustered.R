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
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the integrated data
path_rds <- "./Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/subset_recluster/subset_C3N-00317-T1_immune_unknown_and_recluster/20210305.v1/C3N-00317-T1.Immune_Unknown.Reclustered.20210305.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")

# specify the genes to plot -----------------------------------------------
genes_plot <- c("PTPRC", "CA9", "CD2", "CD3D", "CD3E", "CD3G", "TRAC", "IL7R", "NKG7", "PRF1")
genes_plot <- c("CD8A", "CD8B")
genes_plot <- c("ITGAM", "ITGAX", "CD74")
genes_plot <- c("GZMA")
genes_plot <- c("CD4")

# plot by gene ------------------------------------------------------------
DefaultAssay(srat) <- "RNA"
for (gene_plot in genes_plot) {
  p <- FeaturePlot(object = srat, features = gene_plot, order = T, 
                   cols = c("yellow", "red"))
  ## write output
  file2write <- paste0(dir_out, gene_plot, ".png")
  png(file2write, width = 600, height = 500, res = 150)
  print(p)
  dev.off()
}

for (gene_plot in genes_plot) {
  p <- FeaturePlot(object = srat, features = gene_plot, order = T, 
                   cols = c("yellow", "red"), split.by = "seurat_clusters")
  ## write output
  file2write <- paste0(dir_out, gene_plot, ".split.png")
  png(file2write, width = 3000, height = 400, res = 150)
  print(p)
  dev.off()
}
