# Yige Wu @WashU May 2022
## findmarkers based on the table S2 cluster annotation
## just use cells from the ccRCC patient and only tumor tissue samples (instead of using normal tissue)

# set up working directory and libraries --------------------------------------------------------
## getting the path to the current script - need to keep in the front
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
plan("multiprocess", workers = 8)
options(future.globals.maxSize = 10000 * 1024^2)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)


# input data ------------------------------------------------------
path_rds <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/Analysis_Results/process_external/process_Zhang_scRNA_PNAS_2021/generate_seurat_object_merged_katmai/20220624.v1/Zhang_scRNA_PNAS_2021.Merged.20220624.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading the RDS file!\n")

# set parameters ----------------------------------------------------------
logfc.threshold.run <- 0
min.pct.run <- 0
min.diff.pct.run <- 0
sources_group1 <- c("18854","18855","19703","22368","22604","23459","23843")
sources_group2 <- c("18856","19704","21255","21256","22369","22605")
## group 1
### anno column with "Tumor"
clusterids_group1 <- c("Tumor"); clusterids_group1 
## group 2
clusterids_group2 <- unique(srat@meta.data$anno[srat@meta.data$orig.ident %in% sources_group2]); clusterids_group2 <- clusterids_group2[clusterids_group2 != "unknown"]; clusterids_group2

# process -----------------------------------------------------------------
DefaultAssay(srat)<-"RNA"
dim(srat)
srat@meta.data$cell_group_process <- paste0(srat@meta.data$Source, "_", srat@meta.data$ClusterID)
Idents(srat) <- "cell_group_process"
table(Idents(srat))
idents_group1 <- paste0(sources_group1, "_", clusterids_group1); idents_group1 <- idents_group1[idents_group1 %in% unique(Idents(srat))]; idents_group1
idents_group2 <- paste0(sources_group2, "_", clusterids_group2); idents_group2 <- idents_group2[idents_group2 %in% unique(Idents(srat))]; idents_group2

markers <- FindMarkers(object = srat, test.use = "wilcox", ident.1 = idents_group1, ident.2 = idents_group2, only.pos = F,
                       min.pct = min.pct.run, logfc.threshold = logfc.threshold.run, min.diff.pct = min.diff.pct.run, verbose = T)
markers$gene_symbol <- rownames(markers)

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "logfcthreshold.", logfc.threshold.run, ".minpct.", min.pct.run, ".mindiffpct.", min.diff.pct.run, ".tsv")
write.table(x = markers, file = file2write, quote = F, sep = "\t", row.names = F)


