# Yige Wu @WashU May 2021
## source activate ccrcc_snrna

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
## library additional libaries
library(Signac)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the merged object
atac <- readRDS("./Resources/Analysis_Results/snatac/merge_objects/merge_with_overlapped_peaks_786O_celllines_katmai/20210603.v1/786O_CellLines.PeakOverlapped.20210528.v1.RDS")
## specify the peak to plot
peak_plot <- c("chr9-116234026-116234526")

# process peak ------------------------------------------------------------
chr=strsplit(x = peak_plot, split = "\\-")[[1]][1]
st=strsplit(x = peak_plot, split = "\\-")[[1]][2]; st = as.numeric(st)
en=strsplit(x = peak_plot, split = "\\-")[[1]][3]; en = as.numeric(en)
new_st=st-1000
new_en=en+1000
peak_plot_expanded=paste(chr,new_st,new_en,sep='-')
## change atac ident
print(head(atac@meta.data))
# plot --------------------------------------------------------------------
p=Signac::CoveragePlot(
  object = atac, group.by = "Piece_ID",
  region = peak_plot_expanded,
  annotation = TRUE,
  peaks = F, ranges = Signac::StringToGRanges(peak_plot, sep = c("-", "-")),
  links=FALSE)
## write output
file2write <- paste0(dir_out, "test.png")
png(file2write, width = 1000, height = 1000, res = 150)
print(p)
dev.off()

