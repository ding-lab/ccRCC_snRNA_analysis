# Yige Wu @WashU May 2021
## source activate ccrcc_snrna
## reference: https://satijalab.org/signac/articles/merging.html

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
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
## library additional libaries
library(Signac)
library(Seurat)
library(GenomicRanges)
library(future)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)


# set up running parameters -----------------------------------------------
###some parallelization-solution from the tutorial:
future::plan("multiprocess", workers = 40)
options(future.globals.maxSize = 300 * 1024^3) # for 300 Gb RAM
print("Finished setting running parameters")

# input individual RDS objects into a list ---------------------------------------------------
## input the merged object
atac <- readRDS(file = "./Resources/Analysis_Results/snatac/merge_objects/merge_objects_786O_celllines/20210527.v1/786O_CellLines.Merged.20210527.v1.RDS")
print("Finished readRDS")

# Create a gene activity matrix -------------------------------------------
gene.activities <- GeneActivity(atac)
atac[['RNA']] <- CreateAssayObject(counts = gene.activities)
atac <- NormalizeData(
  object = atac,
  assay = 'RNA',
  normalization.method = 'LogNormalize',
  scale.factor = median(atac$nCount_RNA)
)

## save output
file2write <- paste0(dir_out, "786O_CellLines.Merged.", run_id, ".RDS")
saveRDS(object = atac, file = file2write, compress = T)
print("Finished saveRDS")


