# Yige Wu @WashU Aug 2020

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
## input the integrated data
path_rds <- "./Resources/Analysis_Results/recluster/recluster_cell_groups_in_integrated_data/subset_and_recluster/subset_stroma_and_recluster_on_katmai/20201005.v1/Stroma.Merged.20201005.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")

# fetch data  ----------------------------------------------------
metadata_df <- FetchData(object = srat, vars = c("orig.ident", "UMAP_1", "UMAP_2", "original_barcode"))

## get the genes within the cell type marker table
## save plot
file2write <- paste0(dir_out, "stroma_reclustered", ".metadata.", run_id, ".tsv")
write.table(x = metadata_df, file = file2write, quote = F, sep = "\t", row.names = F)
print("Finish writing the output")


