# Yige Wu @WashU Aug 2020

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
library(tibble)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the integrated data
path_rds <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA//Resources/Analysis_Results/seurat_processing/process_Lake_et_al_katmai/20210819.v1/Lake_etal.Seurat.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")

# find variable genes -----------------------------------------------------
cat("Start running FindVariableFeatures\n")
DefaultAssay(object = srat) <- "SCT"
srat <- FindVariableFeatures(object = srat, selection.method = "vst", nfeatures = 3000, verbose = T)
var_features_df <- data.frame(gene = srat@assays$RNA@var.features)
cat("###########################################\n")

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "findvariablefeatures.", "Lake_etal", ".", run_id, ".tsv")
write.table(var_features_df, file = file2write, quote = F, sep = "\t", row.names = F)
cat("Finished saving the output\n")
cat("###########################################\n")


