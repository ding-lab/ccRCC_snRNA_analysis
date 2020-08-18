# Yige Wu @WashU Feb 2020

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
srat <- readRDS(file = "./Resources/Analysis_Results/integration/merge_different_cases/merge_tumorcells_C3L-00010_C3L-00416_C3L-00583_C3N-00242_on_katmai/20200817.v1/C3L-00416_C3L-00583_C3L-00010_C3N-00242.Tumor_Segments.Merged.20200817.v1.RDS")
print("Finish reading the RDS file!")

# set parameters for findmarkers ------------------------------------------
## set min.pct
logfc.threshold.wilcox <- 0.1
min.pct.wilcox <- 0.1

# specify the two clusters ------------------------------------------------
cluster1 <- 0
cluster2 <- c(2, 4, 5, 6, 8, 9, 10, 11, 12)

# find all markers --------------------------------------------------------
markers_wilcox_df <- FindMarkers(object = srat, test.use = "wilcox", only.pos = F, min.pct = min.pct.wilcox, logfc.threshold = logfc.threshold.wilcox, verbose = T, 
                                 ident.1 = cluster1, ident.2 = cluster2)
markers_wilcox_df$gene <- rownames(markers_wilcox_df)
markers_wilcox_df$ident1 <- cluster1
markers_wilcox_df$ident2 <- paste0(cluster2, collapse = ";")

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "C3N-00010.correlatedtumorcells.FindAllMarkers.newcluster0vsC3L-00416clusters.Wilcox.Minpct", min.pct.wilcox, ".Logfc", logfc.threshold.wilcox,".tsv")
write.table(markers_wilcox_df, file = file2write, 
            quote = F, sep = "\t", row.names = F)