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
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  # "limma",
  "Seurat",
  "future",
  "future.apply"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
# set up future for parallelization
plan("multiprocess", workers = 25)
options(future.globals.maxSize = 10000 * 1024^2)
# registerDoParallel(cores = 6)


# input dependencies ------------------------------------------------------
path_rds <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/Analysis_Results/integration/seuratintegrate_34_ccRCC_samples/reciprocalPCA_integrate_30_ccRCC_tumorcells/20220404.v1/ccRCC.34samples.Tumorcells.SeuratIntegrated.20220404.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading the RDS file!\n")
logfc.threshold.run <- 0
min.pct.run <- 0
min.diff.pct.run <- 0
## input the barcode-to-cluster results
# barcode2cluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/integration/seuratintegrate_34_ccRCC_samples/FindClusters_30_ccRCC_tumorcells_changeresolutions/20220405.v1/ccRCC.34Sample.Tumorcells.Integrated.ReciprocalPCA.Metadata.ByResolution.20220405.v1.tsv")
barcode2cluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/integration/seuratintegrate_34_ccRCC_samples/FindClusters_30_ccRCC_tumorcells_changeresolutions/20220408.v1/ccRCC.34Sample.Tumorcells.Integrated.ReciprocalPCA.Metadata.ByResolution.20220408.v1.tsv")

# process -----------------------------------------------------------------
cat(paste0("running FindMarkers for resolution = 1", "!\n"))
srat@meta.data$cluster_test <- mapvalues(x = rownames(srat@meta.data), from = barcode2cluster_df$barcode, to = as.vector(barcode2cluster_df$integrated_snn_res.1))
Idents(srat) <- "cluster_test"
idents_1 <- c(0, 3, 7, 12)
idents_2 <- unique(barcode2cluster_df$integrated_snn_res.1); idents_2 <- idents_2[!(idents_2 %in% idents_1)]

markers <- FindMarkers(object = srat, test.use = "wilcox", only.pos = F, ident.1 = idents_1, ident.2 = idents_2,
                          min.pct = min.pct.run, logfc.threshold = logfc.threshold.run, min.diff.pct = min.diff.pct.run, verbose = T)
markers$gene_symbol <- rownames(markers)


# write output ------------------------------------------------------------
## set run id
run_id <- "0_3_7_12_vs_others.v2"
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out_parent <- makeOutDir_katmai(path_this_script)
dir_out <- paste0(dir_out_parent, run_id, "/")
dir.create(dir_out)
path_markers <- paste0(dir_out, "res.1", ".tumorcellsreclustered.markers.logfcthreshold.", 
                       logfc.threshold.run, ".minpct.", min.pct.run, ".mindiffpct.", min.diff.pct.run, 
                       ".", run_id, ".tsv")
write.table(x = markers, file = path_markers, quote = F, sep = "\t", row.names = F)

