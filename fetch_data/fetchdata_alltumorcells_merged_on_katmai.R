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
srat <- readRDS(file = "./Resources/Analysis_Results/recluster/recluster_cell_groups_in_integrated_data/subset_and_recluster/subset_tumorcells_and_recluster_on_katmai/20200824.v1/AllTumorCells.Merged.20200824.v1.RDS")
print("Finish reading the RDS file!")

# fetch data --------------------------------------------------------------
umap_data <- Seurat::FetchData(object = srat, vars = c("orig.ident", "ident", "UMAP_1", "UMAP_2"))
umap_data$barcode <- rownames(umap_data)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "AllTumorCells.Merged.Metadata.", run_id, ".tsv")
write.table(x = umap_data, file = file2write, quote = F, sep = "\t", row.names = F)


