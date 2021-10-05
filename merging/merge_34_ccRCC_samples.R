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
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)
options(future.globals.maxSize = 1000 * 1024^2)

# input dependencies ------------------------------------------------------
## input seurat paths
paths_srat <- fread(data.table = F, input = "./Resources/Analysis_Results/data_summary/write_individual_srat_object_paths/20210729.v1/Seurat_Object_Paths.20210729.v1.tsv")
## input doublet prediction
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20210729.v1/scrublet.united_outputs.20210729.v1.tsv", data.table = F)

## scale data with all the features
path_final_file <- paste0(dir_out, "ccRCC.34samples.Merged.", run_id, ".RDS")
path_sct_file <- paste0(dir_out, "SCT." , ".RDS")

if (!file.exists(path_sct_file)) {
  # input the seurat object and store in a list--------------------------------------------------------
  paths_srat2process <- paths_srat %>%
    filter(Case != "C3L-00359")
  renal.list <- list()
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
    renal.list[[sample_id_tmp]] <- seurat_sub_obj
  }
  length(renal.list)
  rm(seurat_obj)
  
  # Run the standard workflow for visualization and clustering ------------
  renal.list <- lapply(X = renal.list, FUN = function(x) {
    x <- NormalizeData(x)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
  })
  ## integrate without anchor
  renal.integrated <- merge(x = renal.list[[1]], y = renal.list[2:length(renal.list)], project = "integrated")
  rm(renal.list)
  
  renal.integrated <- SCTransform(renal.integrated, vars.to.regress = c("nCount_RNA","percent.mito"), return.only.var.genes = F)
  cat("Finished SCTransform!\n")
  saveRDS(object = renal.integrated, file = path_sct_file, compress = T)
  cat("Finished Writing SCTransform!\n")
} else {
  renal.integrated <- readRDS(file = path_sct_file)
}

## keep it consistant with individual processing pipeline
renal.integrated <- RunPCA(renal.integrated, npcs = 30, verbose = FALSE)
cat("Finished RUNPCA!\n")
renal.integrated <- RunUMAP(renal.integrated, reduction = "pca", dims = 1:30)
cat("Finished RUNUMAP!\n")
renal.integrated <- FindNeighbors(renal.integrated, reduction = "pca", dims = 1:30, force.recalc = T)
cat("Finished FindNeighbors!\n")
renal.integrated <- FindClusters(renal.integrated, resolution = 0.5)
cat("Finished FindClusters!\n")
## save as RDS file
saveRDS(object = renal.integrated, file = path_final_file, compress = T)
cat("Finished saving the output!\n")
