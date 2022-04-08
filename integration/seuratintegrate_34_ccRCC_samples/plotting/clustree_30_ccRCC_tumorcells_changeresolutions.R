# Yige Wu @WashU Apr 2022
## https://cran.r-project.org/web/packages/clustree/vignettes/clustree.html

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
  "future.apply",
  "clustree"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
# set up future for parallelization
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 10000 * 1024^2)
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
path_rds <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/Analysis_Results/integration/seuratintegrate_34_ccRCC_samples/reciprocalPCA_integrate_30_ccRCC_tumorcells/20220404.v1/ccRCC.34samples.Tumorcells.SeuratIntegrated.20220404.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading the RDS file!\n")
## input the barcode-to-cluster results
barcode2cluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/integration/seuratintegrate_34_ccRCC_samples/FindClusters_30_ccRCC_tumorcells_changeresolutions/20220405.v1/ccRCC.34Sample.Tumorcells.Integrated.ReciprocalPCA.Metadata.ByResolution.20220405.v1.tsv")

# # FindClusters -----------------------------------------------------------------
# srat <- FindClusters(srat, resolution = c(0.1, 0.2, 0.3, 0.4, 0.5, 1, 2))
# cat("Finished FindClusters!\n")

# add cluster id ----------------------------------------------------------
columns_process <- colnames(barcode2cluster_df)[grepl(pattern = "integrated_snn_res", x = colnames(barcode2cluster_df))]
columns_process <- columns_process[!(columns_process %in% colnames(srat@meta.data))]
columns_process <- c("barcode", columns_process)
meta.data_df <- srat@meta.data
meta.data_df$barcode <- rownames(meta.data_df)
head(meta.data_df)
meta.data_df <- merge(x = meta.data_df, y = barcode2cluster_df[, columns_process], by = c("barcode"), all.x = T, sort = F)
rownames(meta.data_df) <- meta.data_df$barcode
head(meta.data_df)
srat@meta.data <- meta.data_df

# clustree --------------------------------------------------------------
file2write <- paste0(dir_out, "clustree.png")
png(file2write, width = 1500, height = 1300, res = 150)
clustree(srat, prefix = "integrated_snn_res.")
dev.off()
