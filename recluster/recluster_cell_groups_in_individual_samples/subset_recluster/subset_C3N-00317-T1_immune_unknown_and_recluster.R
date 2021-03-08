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
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the integrated data
path_rds <- "./Resources/snRNA_Processed_Data/scRNA_auto/outputs/CPT0012280004/pf3000_fmin200_fmax10000_cmin3000_cmax10000_mito_max0.1/CPT0012280004_processed.rds"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")
## name of the cell group to recluster
table(srat@meta.data$seurat_clusters)
clusters2process <- c(2, 4, 5, 6, 10, 11)
cat(paste0("Name of the cell clusters to recluster: ", clusters2process, "\n"))
cat("###########################################\n")
resolution_tmp <- 2

# subset ------------------------------------------------------------------
## subset
srat <- subset(srat, idents = clusters2process)
dim(srat)

# reclustering ------------------------------------------------------------
## get variably expressed genes
cat("Start running SCTransform\n")
srat <- SCTransform(srat, vars.to.regress = c("nCount_RNA","percent.mito"), return.only.var.genes = F)
cat("###########################################\n")
## keep it consistant with individual processing pipeline
cat("Start running RunPCA\n")
srat <- RunPCA(srat, npcs = num_pc, verbose = FALSE)
cat("###########################################\n")
cat("Start running RunUMAP\n")
srat <- RunUMAP(srat, reduction = "pca", dims = 1:num_pc)
cat("###########################################\n")
cat("Start running FindNeighbors\n")
srat <- FindNeighbors(srat, reduction = "pca", dims = 1:num_pc, force.recalc = T)
cat("###########################################\n")
cat("Start running FindClusters\n")
srat <- FindClusters(srat, resolution = resolution_tmp)
cat("###########################################\n")
## save output
cat("Start saving the reclustered seurat object\n")
file2write <- paste0(dir_out, "C3N-00317-T1.", "Immune_Unknown.", "Reclustered.", "Res", resolution_tmp, ".RDS")
saveRDS(object = srat, file = file2write, compress = T)
cat("###########################################\n")

# write dimplot -----------------------------------------------------------
colors_cluster <- Polychrome::dark.colors(n = length(unique(srat@meta.data$seurat_clusters)))
names(colors_cluster) <- 0:max(as.numeric(as.vector(srat@meta.data$seurat_clusters)))
p <- DimPlot(object = srat)
p <- p + scale_color_manual(values = colors_cluster)
file2write <- paste0(dir_out, "dimplot.png")
png(file2write, width = 800, height = 700, res = 150)
print(p)
dev.off()
