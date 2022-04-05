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
  "future",
  "future.apply"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
# set up future for parallelization
plan("multiprocess", workers = 4)
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
logfc.threshold.run <- 0.25
min.pct.run <- 0.1
min.diff.pct.run <- 0

# process -----------------------------------------------------------------
Idents(srat) <- "integrated_snn_res.0.5"
clusters <- unique(Idents(srat))
pairwise <- combn(clusters, 2)

results_df <- NULL
for (i in 1:ncol(pairwise)) {
  markers <- FindMarkers(object = srat, test.use = "wilcox", ident.1 = pairwise[1, i], ident.2 = pairwise[2, i], only.pos = T,
                         min.pct = min.pct.run, logfc.threshold = logfc.threshold.run, min.diff.pct = min.diff.pct.run, verbose = T)
  markers$gene_symbol <- rownames(markers)
  markers$ident.1 <- pairwise[1, i]
  markers$ident.2 <- pairwise[2, i]
  results_df <- rbind(results_df, markers)
}

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "tumorcellsreclustered.pairwisebycluster.markers.logfcthreshold.", logfc.threshold.run, ".minpct.", min.pct.run, ".mindiffpct.", min.diff.pct.run, ".tsv")
write.table(x = results_df, file = file2write, quote = F, sep = "\t", row.names = F)


