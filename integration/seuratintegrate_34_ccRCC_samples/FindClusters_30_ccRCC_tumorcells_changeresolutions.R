# Yige Wu @WashU Mar 2021
## https://satijalab.org/seurat/archive/v3.0/integration.html
## also used references

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
# dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "Seurat",
  "future",
  "future.apply"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
# set up future for parallelization
plan("multiprocess", workers = 6)
options(future.globals.maxSize = 10000 * 1024^2)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)


# input dependencies ------------------------------------------------------
path_rds <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/Analysis_Results/integration/seuratintegrate_34_ccRCC_samples/reciprocalPCA_integrate_30_ccRCC_tumorcells/20220404.v1/ccRCC.34samples.Tumorcells.SeuratIntegrated.20220404.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading the RDS file!\n")

# process -----------------------------------------------------------------
srat <- FindClusters(srat, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 1, 2, 3, 4))
cat("Finished FindClusters!\n")

# fetch data --------------------------------------------------------------
metadata_df <- srat@meta.data
metadata_df$barcode <- rownames(metadata_df)

file2write <- paste0(dir_out, "ccRCC.34Sample.Tumorcells.Integrated.ReciprocalPCA.Metadata.ByResolution.", run_id, ".tsv")
write.table(x = metadata_df, file = file2write, quote = F, sep = "\t", row.names = F)
print("Finish writing the output!\n")
