# Yige Wu @WashU Sep 2021
## because cluster 10 shows both mast cell and macrophage markers

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
srat <- readRDS(file = "../ccRCC_ST/Processed_Data/Seurat/Outputs/TWFU-HT282N1-S1H3A3N1Z1_1Bmn1_1_5.0/TWFU-HT282N1-S1H3A3N1Z1_1Bmn1_1_5.0_processed_multiomic.rds")
print("Finish reading RDS file")
## name of the cell group to recluster
table(srat@meta.data$seurat_clusters)
clusters2process <- c(10)
cat(paste0("Name of the cell clusters to recluster: ", clusters2process, "\n"))
cat("###########################################\n")
resolution_tmp <- 2
num_pc <- 30

# subset ------------------------------------------------------------------
## subset to non-doublet
aliquot_show <- "HT282N1_S1H3A3N1Z1"
## subset to selected clusters
srat <- subset(srat, idents = clusters2process)
dim(srat)

# reclustering ------------------------------------------------------------
DefaultAssay(srat) <- "RNA"
## get variably expressed genes
cat("Start running SCTransform\n")
srat <- SCTransform(srat, vars.to.regress = c("nCount_RNA","percent.mt","S.Score", "G2M.Score"), return.only.var.genes = F)
cat("###########################################\n")
## keep it consistant with individual processing pipeline
cat("Start running RunPCA\n")
srat <- RunPCA(srat, npcs = num_pc, verbose = FALSE)
cat("###########################################\n")
cat("Start running RunUMAP\n")
srat <- RunUMAP(object = srat, dims = 1:num_pc, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
cat("###########################################\n")
cat("Start running FindNeighbors\n")
srat <- FindNeighbors(srat, reduction = "pca", dims = 1:num_pc, force.recalc = T)
cat("###########################################\n")
cat("Start running FindClusters\n")
srat <- FindClusters(srat, resolution = resolution_tmp)
cat("###########################################\n")
## save output
cat("Start saving the reclustered seurat object\n")
file2write <- paste0(dir_out, aliquot_show, ".C10.", "RNAReclustered.", "Res", resolution_tmp, ".RDS")
saveRDS(object = srat, file = file2write, compress = T)
cat("###########################################\n")

# write dimplot -----------------------------------------------------------
colors_cluster <- Polychrome::dark.colors(n = length(unique(srat@meta.data$seurat_clusters)))
names(colors_cluster) <- 0:max(as.numeric(as.vector(srat@meta.data$seurat_clusters)))
p <- DimPlot(object = srat,  reduction = "umap.rna", group.by = "seurat_clusters", label = TRUE, label.size = 2.5, repel = TRUE)
p <- p + scale_color_manual(values = colors_cluster)
file2write <- paste0(dir_out, "RNA.dimplot.png")
png(file2write, width = 800, height = 700, res = 150)
print(p)
dev.off()
