# Yige Wu @WashU Apr 2021

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
## input seurat paths
paths_srat_df <- fread(data.table = F, input = "./Data_Freezes/V1/snRNA/Tumor_Cell_Reclustered/Paths_TumorCellOnlyReclustered_SeuratObject.20201119.v1.tsv")

# specify genes to plot ---------------------------------------------------
genes_plot <- c("CA9")

easyid_tmp <- "C3N-00733-T1"
path_srat_katmai <- paths_srat_df$Path_katmai[paths_srat_df$Aliquot.snRNA.WU == easyid_tmp]
path_srat_relative <- gsub(x = path_srat_katmai, pattern = "\\/diskmnt\\/Projects\\/ccRCC_scratch\\/ccRCC_snRNA\\/", replacement = "./")
## input seurat object,
srat <- readRDS(file = path_srat_relative)
srat <- readRDS(file = "./Data_Freezes/V1/snRNA/Tumor_Cell_Reclustered/CPT0025880013")

DefaultAssay(srat) <- "RNA"

# plot by gene ------------------------------------------------------------
for (gene_plot in genes_plot) {
  p <- FeaturePlot(object = srat, features = gene_plot, order = T, min.cutoff = "q1", max.cutoff = "q99", cols = c("blue", "red"), label = T)
  ## write output
  file2write <- paste0(dir_out, gene_plot, ".png")
  png(file2write, width = 1000, height = 1000, res = 150)
  print(p)
  dev.off()
}

