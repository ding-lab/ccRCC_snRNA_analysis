# Yige Wu @WashU Mar 2021
## for merging 32 snRNA datasets

# set up libraries and output directory -----------------------------------
## getting the path to the current script
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
path_this_script <- thisFile()
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
# version_tmp <- 1
# run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
run_id <- "20220221.v1"
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)
options(future.globals.maxSize = 1000 * 1024^2)

# input dependencies ------------------------------------------------------
## input seurat paths
paths_srat <- fread(data.table = F, input = "./Resources/Analysis_Results/data_summary/write_individual_srat_object_paths/20210729.v1/Seurat_Object_Paths.20210729.v1.tsv")
## input doublet prediction
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20210729.v1/scrublet.united_outputs.20210729.v1.tsv", data.table = F)

# process -----------------------------------------------------------------
path_anchor_file <- paste0(dir_out, "anchor" , ".RDS")

if (!file.exists(path_anchor_file)) {
  paths_srat2process <- paths_srat %>%
    filter(Case != "C3L-00359")
  srat_list <- list()
  for (i in 1:nrow(paths_srat2process)) {
    sample_id_tmp <- paths_srat2process$Aliquot[i]
    seurat_obj_path <- paths_srat2process$Path_katmai_seurat_object[i]
    seurat_obj <- readRDS(file = seurat_obj_path)
    ## take out the doublets
    barcode2scrublet_tmp_df <- barcode2scrublet_df %>%
      filter(Aliquot == sample_id_tmp) %>%
      filter(Barcode %in% rownames(seurat_obj@meta.data)) %>%
      filter(!predicted_doublet)
    barcodes_keep <- barcode2scrublet_tmp_df$Barcode
    ###
    print("subsetting")
    seurat_sub_obj <- subset(x = seurat_obj, cells = barcodes_keep)
    print(dim(seurat_sub_obj))
    srat_list[[sample_id_tmp]] <- seurat_sub_obj
  }
  length(srat_list)
  rm(seurat_obj)
  cat("Finished making seurat object list!\n")
  
  # Run the standard workflow for visualization and clustering ------------
  srat_list <- lapply(X = srat_list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
  })
  cat("Finished NormalizeData and FindVariableFeatures for the list!\n")
  
  # select features that are repeatedly variable across datasets for integration
  features <- SelectIntegrationFeatures(object.list = srat_list)
  cat("Finished SelectIntegrationFeatures!\n")
  
  srat_anchors <- FindIntegrationAnchors(object.list = srat_list, anchor.features = features)
  cat("Finished FindIntegrationAnchors!\n")
  
  saveRDS(object = srat_anchors, file = path_anchor_file, compress = T)
  cat("Finished saveRDS srat_anchors!\n")
}  else {
  srat_anchors <- readRDS(file = path_anchor_file)
}

# this command creates an 'integrated' data assay
srat_integrated <- IntegrateData(anchorset = srat_anchors)
cat("Finished IntegrateData!\n")
rm(srat_list)
rm(srat_anchors)
cat("Finished deleting list and anchor!\n")

## keep it consistant with individual processing pipeline
DefaultAssay(srat_integrated) <- "integrated"
srat_integrated <- ScaleData(srat_integrated, verbose = F)
cat("Finished ScaleData!\n")
# srat_integrated <- RunPCA(srat_integrated, npcs = 30, verbose = FALSE)
srat_integrated <- RunPCA(srat_integrated, npcs = 45, verbose = FALSE)
cat("Finished RUNPCA!\n")
srat_integrated <- RunUMAP(srat_integrated, reduction = "pca", dims = 1:30)
cat("Finished RUNUMAP!\n")
srat_integrated <- FindNeighbors(srat_integrated, reduction = "pca", dims = 1:30, force.recalc = T)
cat("Finished FindNeighbors!\n")
srat_integrated <- FindClusters(srat_integrated, resolution = 0.5)
cat("Finished FindClusters!\n")
## save as RDS file
## scale data with all the features
path_final_file <- paste0(dir_out, "ccRCC.34samples.SeuratIntegrated.", run_id, ".RDS")
saveRDS(object = srat_integrated, file = path_final_file, compress = T)
cat("Finished saving the output!\n")
