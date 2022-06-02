# Yige Wu @WashU Jun 2022
## finding differentially expressed gene for each cell type using integrared object

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
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "Seurat",
  "Signac",
  "ggplot2",
  "ggpubr"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    print(paste0(pkg_name_tmp, "is being installed!"))
    BiocManager::install(pkgs = pkg_name_tmp, update = F)
    install.packages(pkg_name_tmp, dependencies = T)
  }
  print(paste0(pkg_name_tmp, " is installed!"))
  library(package = pkg_name_tmp, character.only = T)
}
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input data ------------------------------------------------------
## input the integrated data
path_rds <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/Signac.1.0.0/3.Merge_snATAC/Merge.SelectPeaks.v.20210706/28_ccRCC_snATAC.selectedPeaks.chromvar.cicero.v3.20210725.rds"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")
## input the barcode-cell-type table
barcode2celltype_df <- fread(input = "./Resources/snATAC_Processed_Data/Barcode_Annotation/28_ccRCC_snATAC_ManualReviwed.v2.20210709.tsv", data.table = F)
cat("finish reading the barcode-to-cell type table!\n")

# specify the parameters --------------------------------------------------
motif_id_plot <- "MA1107.2"

# make plot data ----------------------------------------------------
# this will extract chromvar assay:
DefaultAssay(srat) <- 'chromvar'
chromv = GetAssayData(object = srat)
## make plot data
chromv_df <- data.frame(barcode_uniq = colnames(chromv), motif_score = as.vector(chromv[motif_id_plot,]))
chromv_df$cell_group_plot <- mapvalues(x = chromv_df$barcode_uniq , from = barcode2celltype_df$V1, to = as.vector(barcode2celltype_df$cell_type))
plotdata_df <- chromv_df %>%
  filter(cell_group_plot %in% c("Tumor", "EMT tumor cells", "PT"))
## write output
file2write <- paste0(dir_out, motif_id_plot, ".motif_score_by_cell.tsv")
write.table(x = plotdata_df, file = file2write, quote = F, sep = "\t", row.names = F)




