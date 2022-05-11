# Yige Wu @WashU May 2022

# set up working directory and libraries --------------------------------------------------------
## set working directory
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
setwd(dir_base)
## load libraries
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
# set up future for parallization
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 10000 * 1024^2)

# set up output directory -----------------------------------
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
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
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input data ------------------------------------------------------
path_rds <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/Analysis_Results/process_external/process_Young_scRNA_Science_2018/generate_seurat_object_merged_katmai/20220510.v1/Young_scRNA_Science_2018.Merged.20220510.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading the RDS file!\n")

# set parameters ----------------------------------------------------------
logfc.threshold.run <- 0
min.pct.run <- 0
min.diff.pct.run <- 0
idents_group1 <- c("T4", "T6", "T7", "T9", "T10", "T12", "T17"); idents_group1
idents_group2 <- c(paste0("IT", c(0, 2, 3, 4, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24)),
                   "T0", "T1", "T2", "T3", "T8", "T11", "T14"); idents_group2

# process -----------------------------------------------------------------
Idents(srat) <- "ClusterID"
table(Idents(srat))
markers <- FindMarkers(object = srat, test.use = "wilcox", ident.1 = idents_group1, ident.2 = idents_group2, only.pos = F,
                       min.pct = min.pct.run, logfc.threshold = logfc.threshold.run, min.diff.pct = min.diff.pct.run, verbose = T)
markers$gene_symbol <- rownames(markers)

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "logfcthreshold.", logfc.threshold.run, ".minpct.", min.pct.run, ".mindiffpct.", min.diff.pct.run, ".tsv")
write.table(x = markers, file = file2write, quote = F, sep = "\t", row.names = F)


