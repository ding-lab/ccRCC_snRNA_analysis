# Yige Wu @WashU Nov 2020

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
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(Seurat)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input seurat object paths
srat_paths_df <- fread(data.table = F, input = "./Data_Freezes/V2/snRNA/Tumor_Cell_Reclustered/Paths_TumorCellOnlyReclustered_SeuratObject.20210805.v1.tsv")

# process by each aliquot ----------------------------------------------------
barcode_metadata_df <- NULL
for (aliquot_tmp in srat_paths_df$Aliquot.snRNA.WU) {
  for (gene_exp in c("MYC", "QKI", "ARID1B", "CD70")) {
    ## input srat object
    srat_path <- srat_paths_df$Path_katmai[srat_paths_df$Aliquot.snRNA.WU == aliquot_tmp]
    srat <- readRDS(file = srat_path)
    
    p = FeaturePlot(object = srat, features = gene_exp, max.cutoff = "q90")
    png(paste0(dir_out, aliquot_tmp, ".", gene_exp, ".png"), width = 500, height = 500, res = 150)
    print(p)
    dev.off()
  }
}
print("Finished!")

#