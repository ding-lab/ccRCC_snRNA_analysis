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
# dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
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
srat_paths <- fread(data.table = F, input = "./Resources/Analysis_Results/data_summary/write_individual_srat_object_paths/20210729.v1/Seurat_Object_Paths.20210729.v1.tsv")
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210809.v1/meta_data.20210809.v1.tsv")

# plot --------------------------------------------------------------
aliquots2process <- unique(srat_paths$Aliquot)
pct_df <- NULL

for (aliquot_tmp in aliquots2process) {
  easyid <- idmetadata_df$Aliquot.snRNA.WU[idmetadata_df$Aliquot.snRNA == aliquot_tmp]
  ## input individually processed seurat object
  seurat_obj_path <- srat_paths$Path_katmai_seurat_object[srat_paths$Aliquot == aliquot_tmp]
  seurat_obj_path
  srat <- readRDS(file = seurat_obj_path)
  print(dim(srat))
  
  p <- ElbowPlot(srat, ndims = 40)
  # file2write <- paste0(dir_out, easyid, ".ElbowPlot", ".pdf")
  # pdf(file2write, width = 5, height = 4, useDingbats = F)
  file2write <- paste0(dir_out, easyid, ".ElbowPlot", ".png")
  png(file2write, width = 600, height = 500, res = 150)
  print(p)
  dev.off()
  
  # Determine percent of variation associated with each PC
  pct <- srat[["pca"]]@stdev / sum(srat[["pca"]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  
  pct_tmp_df <- data.frame(easy_id = rep(easyid, length(pct)), rank_pc = 1:length(pct), pct = pct, cumu_pct = cumu, aliquot = rep(aliquot_tmp, length(pct)))
  pct_df <- rbind(pct_tmp_df, pct_df)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  print(co1)
}

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "pct_standard_variances_by_PC.", run_id, ".tsv")
write.table(x = pct_df, file = file2write, quote = F, sep = "\t", row.names = F)


