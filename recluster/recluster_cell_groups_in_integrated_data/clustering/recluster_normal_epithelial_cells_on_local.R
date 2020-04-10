# Yige Wu @WashU Apr 2020
## for each individual sample, isolating tumor cells assigned from the integrated data and re-do clustering

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
srat <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/subset/subset_normal_epithelial_cells/20200406.v1/30_aliquot_integration.20200212.v3.RDS.normal_epithelial_cells.RDS")

# run reclustering by each aliquot ----------------------------------------
## Run the standard workflow for clustering and visualization
srat <- ScaleData(srat, features = rownames(srat@assays$RNA@counts))
srat <- FindVariableFeatures(object = srat, selection.method = "vst", nfeatures = num_features)

## RunPCA
srat <- RunPCA(srat, npcs = num_pc, verbose = FALSE)
srat <- RunUMAP(srat, reduction = "pca", dims = 1:num_pc)
srat <- FindNeighbors(srat, reduction = "pca", dims = 1:num_pc, force.recalc = T)
srat <- FindClusters(srat, resolution = 0.5)

# make table for paths to these objects -----------------------------------
file2write <- paste0(dir_out, "normal_epithelial_cells.reclustered.", run_id, ".RDS")
saveRDS(object = srat, file = file2write, compress = T)