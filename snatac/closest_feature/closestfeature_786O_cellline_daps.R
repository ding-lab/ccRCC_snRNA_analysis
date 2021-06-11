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
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input individual RDS objects into a list ---------------------------------------------------
## input the merged object
atac <- readRDS(file = "./Resources/Analysis_Results/snatac/merge_objects/add_gene_activity_to_786O_celllines_merged_katmai/20210528.v1/786O_CellLines.Merged.20210528.v1.RDS")
print("Finished readRDS")

# query -------------------------------------------------------------------
regions_process <- head(rownames(atac))
print(regions_process)
closest_genes <- ClosestFeature(object = atac, regions = regions_process)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "786O.DA_peaks.Closest_Features.", run_id, ".tsv")
write.table(x = closest_genes, file = file2write, quote = F, row.names = F, sep = "\t")

