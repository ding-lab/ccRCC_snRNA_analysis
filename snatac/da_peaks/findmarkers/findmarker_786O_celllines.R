# Yige Wu @WashU Jun 2021
## source activate ccrcc_snrna

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
library(future)
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)
plan("multiprocess", workers = 10)
options(future.globals.maxSize= 891289600)

# input dependencies ------------------------------------------------------
atac <- readRDS("./Resources/Analysis_Results/snatac/merge_objects/merge_with_overlapped_peaks_786O_celllines_katmai/20210603.v1/786O_CellLines.PeakOverlapped.20210528.v1.RDS")
print(paste0("Finished reading atac object!"))

# specify the run parameters ----------------------------------------------
logfc.threshold.run <- 0
min.pct.run <- 0.1
min.diff.pct.run <- 0

# findmarker --------------------------------------------------------------
print(head(atac@meta.data))
Idents(atac) <- "Piece_ID"
da_peaks <- FindMarkers(
  object = atac,
  ident.1 ='BAP1_786O',
  ident.2='Control_786O',
  only.pos = FALSE,
  min.pct = min.pct.run, min.diff.pct = min.diff.pct.run, logfc.threshold = logfc.threshold.run,
  test.use = 'LR', verbose = T,
  latent.vars = 'peak_RF_500MACS2'
)
da_peaks$row_name <- rownames(da_peaks)
print(paste0("Finished FindMarkers!"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "DA_Peaks.BAP1_vs_Control.786O_CellLines.", 
                     "min.diff.pct.", min.diff.pct.run,
                     "min.pct.", min.pct.run,
                     "logfc.threshold.", logfc.threshold.run, ".tsv")
write.table(x = da_peaks, file = file2write, sep = "\t", row.names = F, quote = F)
print(paste0("Finished write.table!"))


