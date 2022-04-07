require(devtools)
devtools::install_github("YosefLab/VISION")

library(Seurat)
library(VISION)
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
setwd(dir_base)
seurat_object <- readRDS(file = "./Resources/snRNA_Processed_Data/scRNA_auto/outputs/CPT0000870003/pf4000_fmin200_fmax10000_cmin4000_cmax10000_mito_max0.1/CPT0000870003_processed.rds")
signatures <- c("./Resources/Knowledge/Databases/MSigDB/msigdb_v7.4_GMTs/h.all.v7.4.symbols.gmt",
                "./Resources/Knowledge/Databases/MSigDB/msigdb_v7.4_GMTs/c2.cp.v7.4.symbols.gmt",
                "./Resources/Knowledge/Databases/MSigDB/msigdb_v7.4_GMTs/c5.go.bp.v7.4.symbols.gmt")
vision.obj <- Vision(seurat_object, signatures = signatures)

seurat_object <- readRDS(file = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/CPT0000870003/pf4000_fmin200_fmax10000_cmin4000_cmax10000_mito_max0.1/CPT0000870003_processed.rds")
