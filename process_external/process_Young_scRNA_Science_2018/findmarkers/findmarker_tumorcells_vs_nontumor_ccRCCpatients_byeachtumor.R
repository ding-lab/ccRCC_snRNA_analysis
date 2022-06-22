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
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 10000 * 1024^2)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
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
## group 1
### Cell_type1 column with "renal cell carcinoma" value, mostly from RCC1, RCC2, and VHL_RCC, which are the 3 clear cell RCCs
idents_group1 <- c("T4", "T6", "T7", "T9", "T10", "T12", "T17"); idents_group1 
## group 2
##@ T0, 1, 8, 14 are labeled "Normal_cell" in Cell_type1 column, all of them from kidney tumor samples, mostly from RCC1, RCC2, and VHL_RCC, which are the 3 clear cell RCCs
##@ T2, 3, 11 are labeled "Endothelium" in Cell_type1 column, all of them from kidney tumor samples, mostly from RCC1, RCC2, and VHL_RCC, which are the 3 clear cell RCCs
##@ the following IT0-24 are immune cells from kidney tumors, excluding those labeled "junk" in Cell_type1 column
#### except IT4, IT15, IT16, 19 is mostly from papillary; IT17, IT22 mostly from a wilms tumor
idents_group2 <- c(paste0("IT", c(0, 2, 3, 6, 7, 8, 9, 10, 11, 12, 13, 14, 18, 20, 21, 23, 24)),
                   "T0", "T1", "T2", "T3", "T8", "T11", "T14"); idents_group2
sources_process <- c("VHL_Kid_T_ldc_1_1", "VHL_Kid_T_ldc_1_2", "RCC2_Kid_T_ldc_1_1", "RCC2_Kid_T_ldc_1_2", "RCC2_Kid_T_ldc_2_1", "RCC2_Kid_T_ldc_2_2", "RCC1_Kid_T_ldc_1_1", "RCC1_Kid_T_ldc_1_2", "RCC1_Kid_T_ldc_2_1", "RCC1_Kid_T_ldc_2_2")

# process -----------------------------------------------------------------
DefaultAssay(srat)<-"RNA"
dim(srat)
srat@meta.data$cell_group_process <- paste0(srat@meta.data$Source, "_", srat@meta.data$ClusterID)
Idents(srat) <- "cell_group_process"
table(Idents(srat))

# source_tmp <- "VHL_Kid_T_ldc_1_1"
markers <- NULL
for (source_tmp in sources_process) {
  markers_tmp <- FindMarkers(object = srat, test.use = "wilcox", ident.1 = paste0(source_tmp, "_", idents_group1), ident.2 = paste0(source_tmp, "_", idents_group2), only.pos = F,
                         min.pct = min.pct.run, logfc.threshold = logfc.threshold.run, min.diff.pct = min.diff.pct.run, verbose = T)
  markers_tmp$gene_symbol <- rownames(markers_tmp)
  markers_tmp$Source <- source_tmp
  markers <- rbind(markers, markers_tmp)
}

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "logfcthreshold.", logfc.threshold.run, ".minpct.", min.pct.run, ".mindiffpct.", min.diff.pct.run, ".tsv")
write.table(x = markers, file = file2write, quote = F, sep = "\t", row.names = F)


