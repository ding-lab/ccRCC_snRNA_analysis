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
path_rds <- "./Resources/snRNA_Processed_Data/scRNA_auto/outputs/CPT0001260013/pf1000_fmin200_fmax10000_cmin1000_cmax10000_mito_max0.1/CPT0001260013_processed.rds"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")
## input the barcode-manualsubcluster info
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_corrected_by_individual_sample_inspection/20200904.v1/31Aliquot.Barcode2CellType.20200904.v1.tsv", data.table = F)
## name of the cell group to recluster
cellgroup2process <- c("Tumor cells", "Tumor-like cells")
cat(paste0("Name of the cell group to recluster: ", cellgroup2process, "\n"))
cat("###########################################\n")

# subset ------------------------------------------------------------------
## filter barcode2celltype tabl
barcode2celltype_filtered_df <- barcode2celltype_df %>%
  filter(orig.ident == "CPT0001260013") %>%
  filter(Cell_group.detailed %in% cellgroup2process)
nrow(barcode2celltype_filtered_df)
## subset
srat <- subset(srat, cells = barcode2celltype_filtered_df$individual_barcode)
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
