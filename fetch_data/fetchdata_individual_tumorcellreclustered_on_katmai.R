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
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input seurat object paths
srat_paths_df <- fread(data.table = F, input = "./Data_Freezes/V1/snRNA/Tumor_Cell_Reclustered/Paths_TumorCellOnlyReclustered_SeuratObject.20201119.v1.tsv")

# process by each aliquot ----------------------------------------------------
barcode_metadata_df <- NULL
for (aliquot_tmp in srat_paths_df$Aliquot.snRNA.WU) {
  ## input srat object
  srat_path <- srat_paths_df$Path_katmai[srat_paths_df$Aliquot.snRNA.WU == aliquot_tmp]
  srat <- readRDS(file = srat_path)
  
  ## extract current meta data
  barcode_metadata_tmp <- FetchData(object = srat, vars = c("UMAP_1", "UMAP_2", "orig.ident", "seurat_clusters"))
  barcode_metadata_tmp$barcode_tumorcellreclustered <- rownames(barcode_metadata_tmp)
  barcode_metadata_tmp$easy_id <- aliquot_tmp
  
  ## bind with the super table
  barcode_metadata_df <- rbind(barcode_metadata_tmp, barcode_metadata_df)
}

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "MetaData_TumorCellOnlyReclustered.", run_id, ".tsv")
write.table(x = barcode_metadata_df, file = file2write, sep = '\t', quote = F, row.names = F)

  
  