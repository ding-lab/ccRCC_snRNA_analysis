# Yige Wu @WashU Apr 2021

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
## input the integrated data
path_rds <- "/diskmnt/Projects/HTAN_analysis_2/BRCA/Analyses/Dan/Reannotation/HTAN_BRCA_scRNA_v10.rds"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")

## spcify assay
assay_process <- "SCT"
slot_process <- "data"
cat(paste0("Assay: ", assay_process, "\n"))

# set ident ---------------------------------------------------------------
srat@meta.data$Cell_group <- paste0(srat@meta.data$Piece_ID, "_", srat@meta.data$cell_type)
unique(srat@meta.data$Cell_group)
Idents(srat) <- "Cell_group" 

# run average expression --------------------------------------------------
aliquot.averages <- AverageExpression(srat, assays = assay_process, slot = slot_process)
print("Finish running AverageExpression!\n")
cat("###########################################\n")

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "BRCA.", "avgexp.", assay_process, ".", slot_process, ".", "cell_type.byPiece_ID.", run_id, ".tsv")
write.table(aliquot.averages, file = file2write, quote = F, sep = "\t", row.names = T)
cat("Finished saving the output\n")
cat("###########################################\n")


