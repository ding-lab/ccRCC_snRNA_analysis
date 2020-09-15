# Yige Wu @WashU Sep 2020

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
# dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(Signac)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the seurat object
srat <- readRDS(file = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/4.Cell_type_markers/1.ChromVar/out/9_ccRCC_MergedObj_chromVAR.rds")
print("Finish reading the RDS file!")
## input the list of peaks
peaks_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/unite_cellgroup4_da_peaks/20200915.v1/DA_peaks.chromvar.MergedObj.byCell_group4.20200915.v1.tsv")

# get genomic regions to qeury --------------------------------------------
regions_process <- unique(peaks_df$V1)
length(regions_process)

# query -------------------------------------------------------------------
closest_genes <- ClosestFeature(object = srat, regions = regions_process)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "DA_peaks.Closest_Features.", run_id, ".tsv")
write.table(x = closest_genes, file = file2write, quote = F, row.names = F, sep = "\t")


