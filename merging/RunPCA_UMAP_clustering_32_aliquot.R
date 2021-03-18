# Yige Wu @WashU Jul 2020
## for merging 31 snRNA datasets

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
renal.integrated <- readRDS(file = "./Resources/Analysis_Results/merging/SCTransform_32_aliquot/20210318.v1/32_aliquot.SCTransform.20210318.v1.RDS")

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
saveRDS(object = renal.integrated, file = paste0(dir_out, "32_aliquot.Merged.", run_id, ".RDS"), compress = T)
cat("Finished saving the output!\n")
