# Yige Wu @WashU Aug 2020
## for integrating  snRNA datasets belong to the same patient

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
dir_tmp <- getwd()
print(dir_tmp)
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
path_rds <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/Analysis_Results/merging/merge_34_ccRCC_samples/20220218.v1/ccRCC.34samples.Merged.20220218.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading the RDS file!\n")

# plot --------------------------------------------------------------
p <- ElbowPlot(srat, ndims = 45)
file2write <- paste0(dir_out, "ccRCC.34Sample.Merged.ElbowPlot.", run_id, ".pdf")
pdf(file2write, width = 5, height = 4, useDingbats = F)
print(p)
dev.off()

file2write <- paste0(dir_out, "ccRCC.34Sample.Merged.pctvariance.", run_id, ".txt")
sink(file2write)
# Determine percent of variation associated with each PC
pct <- srat[["pca"]]@stdev / sum(srat[["pca"]]@stdev) * 100
print(pct)

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)
print(cumu)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]
print(co1)
sink()
