# Yige Wu @WashU Sep 2020

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
## input the integrated data
path_rds <- "./Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/subset_C3N-01200_tumorlikecells_and_recluster/20200910.v1/TumorLikeCells.Reclustered.20200910.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")

# Dimplot -----------------------------------------------------------------
## make distinguishable colors
head(srat@meta.data$seurat_clusters)
as.numeric(head(srat@meta.data$seurat_clusters))
as.numeric(as.vector(head(srat@meta.data$seurat_clusters)))
colors_cluster <- Polychrome::dark.colors(n = length(unique(srat@meta.data$seurat_clusters)))
names(colors_cluster) <- 0:max(as.numeric(as.vector(srat@meta.data$seurat_clusters)))
p <- DimPlot(object = srat)
p <- p + scale_color_manual(values = colors_cluster)
p
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Seurat_Cluster_Id.", "C3N-01200", ".png")
png(file2write, width = 900, height = 800, res = 150)
print(p)
dev.off()

