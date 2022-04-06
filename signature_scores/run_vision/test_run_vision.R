require(devtools)
install_github("YosefLab/VISION")

library(Seurat)
seurat_object <- readRDS(file = "./Resources/snRNA_Processed_Data/scRNA_auto/outputs/CPT0000870003/pf4000_fmin200_fmax10000_cmin4000_cmax10000_mito_max0.1/CPT0000870003_processed.rds")
vision.obj <- Vision(seurat_object, signatures = signatures)
