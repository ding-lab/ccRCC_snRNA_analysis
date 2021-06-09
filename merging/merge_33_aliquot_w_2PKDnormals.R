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
library(future)
plan("multiprocess", workers = 3)
options(future.globals.maxSize = 35 * 1024^3) # for 7 Gb RAM

# input dependencies ------------------------------------------------------
## input seurat paths
paths_srat <- fread(data.table = F, input = "./Data_Freezes/V2/snRNA/Seurat_Object_Paths.20210428.v1.tsv")
## input doublet prediction
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20210428.v1/scrublet.united_outputs.20210428.v1.tsv", data.table = F)

# make united paths -------------------------------------------------------
paths_srat2process <- paths_srat %>%
  select(Aliquot, Path_katmai_seurat_object)
paths_srat2process <- rbind(data.frame(Aliquot = c("HTHK1900070", 
                                                   "HTHK190043"),
                                       Path_katmai_seurat_object = c("/diskmnt/Projects/Users/rliu/Projects/PKD/Seurat/HTHK1900070/nfeature_400_3000_ncount_1000_7000_mito_0.1/HTHK1900070_processed.rds",
                                                                     "/diskmnt/Projects/Users/rliu/Projects/PKD/Seurat/HTHK190043/nfeature_400_6000_ncount_1000_Inf_mito_0.1/HTHK190043_processed.rds")), 
                            paths_srat2process)

# process -----------------------------------------------------------------
path_final_file <- paste0(dir_out, "33_aliquot_w_2PKDnormals.Merged.", run_id, ".RDS")
path_sct_file <- paste0(dir_out, "33_aliquot_w_2PKDnormals.SCT.", run_id , ".RDS")

if (!file.exists(path_sct_file)) {
  # input the seurat object and store in a list--------------------------------------------------------
  renal.list <- list()
  for (i in 1:nrow(paths_srat2process)) {
    sample_id_tmp <- paths_srat2process$Aliquot[i]
    seurat_obj_path <- paths_srat2process$Path_katmai_seurat_object[i]
    seurat_obj <- readRDS(file = seurat_obj_path)
    if (sample_id_tmp %in% barcode2scrublet_df$Aliquot) {
      ## take out the doublets
      barcode2scrublet_tmp_df <- barcode2scrublet_df %>%
        filter(Aliquot == sample_id_tmp) %>%
        filter(Barcode %in% rownames(seurat_obj@meta.data)) %>%
        filter(!predicted_doublet)
      barcodes_keep <- barcode2scrublet_tmp_df$Barcode
      ### subsetting
      print("subsetting")
      seurat_sub_obj <- subset(x = seurat_obj, cells = barcodes_keep)
      print(dim(seurat_sub_obj))
      renal.list[[sample_id_tmp]] <- seurat_sub_obj
    } else {
      print(dim(seurat_obj))
      renal.list[[sample_id_tmp]] <- seurat_obj
    }

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
