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
path_rds <- "./Resources/snRNA_Processed_Data/Merged_Seurat_Objects/C3N-01200.Tumor_Segments.Merged.20200319.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")
## input the barcode-manualsubcluster info
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_corrected_by_individual_sample_inspection/20200904.v1/31Aliquot.Barcode2CellType.20200904.v1.tsv", data.table = F)
## name of the cell group to recluster
cellgroup2process <- c("Tumor cells", "Tumor-like cells")
cat(paste0("Name of the cell group to recluster: ", cellgroup2process, "\n"))
cat("###########################################\n")

# add cell type to the Seurat meta data---------------------------------------------
BC <- srat@meta.data %>% rownames
## get original barcode
srat@meta.data$original_barcode <- BC %>% strsplit("_") %>% lapply("[[",1) %>% unlist
head(srat@meta.data$original_barcode)
## make combined id for the seurat meta data
srat@meta.data$id_aliquot_barcode <- paste0(srat@meta.data$orig.ident, "_", srat@meta.data$original_barcode)
head(srat@meta.data$id_aliquot_barcode)
## make combined id for the barcode2celltype table
barcode2celltype_df$id_aliquot_barcode <- paste0(barcode2celltype_df$orig.ident, "_", barcode2celltype_df$individual_barcode)
## map cell type shorter
srat@meta.data$Cell_type.shorter <- mapvalues(x = srat@meta.data$id_aliquot_barcode, from = barcode2celltype_df$id_aliquot_barcode, to = as.vector(barcode2celltype_df$Cell_type.shorter))
head(srat@meta.data$Cell_type.shorter)
unique(srat@meta.data$Cell_type.shorter)
Idents(srat) <- "Cell_type.shorter"

# subset ------------------------------------------------------------------
## subset
srat <- subset(srat, idents = cellgroup2process)
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
srat <- FindClusters(srat, resolution = resolution_findclusters)
cat("###########################################\n")
## save output
cat("Start saving the reclustered seurat object\n")
file2write <- paste0(dir_out, "TumorLikeCells.Reclustered.", run_id, ".RDS")
saveRDS(object = srat, file = file2write, compress = T)
cat("###########################################\n")
