# Yige Wu @WashU Mar 2021
## for generating the barcode annotation files (assigned to one random group) for 10Xmapping
## only use barcodes of after QC barcodes so that 10Xmapping won't run forever

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
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
path_cellranger <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/Cell_Ranger/outputs/"

# set parameters ----------------------------------------------------------
aliquots2process <- c("CPT0001180011")

# input seurat object -----------------------------------------------------
for (aliquot in aliquots2process) {
  path_barcodes <- paste0(path_cellranger, aliquot, "/outs/filtered_peak_bc_matrix/barcodes.tsv")
  barcodes_df = read.delim(path_barcodes, header = FALSE,stringsAsFactors = FALSE)
  
  anno_tab_tmp <- barcodes_df %>%
    mutate(random_group = "0")
  write.table(x = anno_tab_tmp, file = paste0(dir_out, aliquot, ".snATAC.", "CellRangerFiltered.Annotation.tsv"), quote = F, row.names = F, sep = "\t", col.names = F)
}


