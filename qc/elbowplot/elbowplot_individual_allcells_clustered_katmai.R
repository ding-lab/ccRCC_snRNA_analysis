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
## input seurat object paths
srat_paths <- fread(data.table = F, input = "./Data_Freezes/V2/snRNA/Seurat_Object_Paths.20210428.v1.tsv")

# plot --------------------------------------------------------------
aliquots2process <- unique(srat_paths$Aliquot)

for (aliquot_tmp in aliquots2process) {
  easyid <- idmetadata_df$Aliquot.snRNA.WU[idmetadata_df$Aliquot.snRNA == aliquot_tmp]
  ## input individually processed seurat object
  seurat_obj_path <- srat_paths$Path_katmai_seurat_object[srat_paths$Aliquot == aliquot_tmp]
  seurat_obj_path
  srat <- readRDS(file = seurat_obj_path)
  print(dim(srat))
  
  p <- ElbowPlot(srat, ndims = 35)
  file2write <- paste0(dir_out, easyid, ".ElbowPlot", ".pdf")
  pdf(file2write, width = 5, height = 5, useDingbats = F)
  print(p)
  dev.off()
  
}



