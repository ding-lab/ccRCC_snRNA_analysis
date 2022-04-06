# Yige Wu @WashU Mar 2022
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
  "VISION"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input  ------------------------------------------------------
## input seurat object
path_rds <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/Analysis_Results/integration/seuratintegrate_34_ccRCC_samples/reciprocalPCA_integrate_34_ccRCC_samples/20220222.v1/ccRCC.34samples.SeuratIntegrated.20220222.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading the RDS file!\n")
## set paths to signature objects
signatures <- c("./Resources/Knowledge/Databases/MSigDB/msigdb_v7.4_GMTs/h.all.v7.4.symbols.gmt",
                "./Resources/Knowledge/Databases/MSigDB/msigdb_v7.4_GMTs/c2.cp.v7.4.symbols.gmt",
                "./Resources/Knowledge/Databases/MSigDB/msigdb_v7.4_GMTs/c5.go.bp.v7.4.symbols.gmt")

# process -----------------------------------------------------------------
## this is important because the default assay for the seurat object is "integrated" and I think if not set to RNA assay Vision will just take integrated assay,
## which will give an error about  "rownames(data) = NULL. Expression matrix must have gene names as the rownames" because
### > rownames(srat$integrated@counts)
### NULL
### rownames(srat$integrated@data will give top variably expressed genes
DefaultAssay(srat) <- "RNA"
vision.obj <- Vision(srat, signatures = signatures)
print("Finish creating the vision object!\n")
# Set the number of threads when running parallel computations
options(mc.cores = 4)
vision.obj <- analyze(vision.obj)
print("Finish analyze the vision object!\n")
sigScores <- getSignatureScores(vision.obj)
print("Finish getSignatureScores!\n")

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "ccRCC.34samples.SeuratIntegrated.Vision.", run_id, ".RDS")
saveRDS(object = vision.obj, file = file2write, compress = T)
file2write <- paste0(dir_out, "ccRCC.34samples.SeuratIntegrated.Vision.scores.", run_id, ".RDS")
saveRDS(object = sigScores, file = file2write, compress = T)
