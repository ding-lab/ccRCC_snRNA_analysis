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
atac <- readRDS(file = "./Resources/Analysis_Results/snatac/merge_objects/merge_objects_786O_celllines/20210527.v1/786O_CellLines.Merged.20210527.v1.RDS")
## specify the peak to plot
peaks_plot_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/BAP1_Specific/BAP1_down_Promoter_and_overlapping_methylationProbes_m1.tsv")
## input the probe info
probes_anno_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/methylation/annotate_methyl_subtype_specific_top1000_probes/20210525.v1/methyl_subtype_specific_1000_probe2gene.20210525.v1.tsv")

# for (i in 1:nrow(peaks_plot_df)) {
for (i in 1) {
    
  # process peak ------------------------------------------------------------
  chr=peaks_plot_df$seqnames[i]
  st=peaks_plot_df$start[i]; st = as.numeric(st)
  en=peaks_plot_df$end[i]; en = as.numeric(en)
  peak_plot=paste(chr,st,en,sep='-')
  new_st=st-1000
  new_en=en+1000
  peak_plot_expanded=paste(chr,new_st,new_en,sep='-')
  gene_tmp <- peaks_plot_df$SYMBOL[i]
  
  # plot --------------------------------------------------------------------
  cov_plot=CoveragePlot(
    object = atac, group.by = "Piece_ID", 
    region = peak_plot_expanded,
    annotation = F,
    peaks = F,
    links=FALSE)
  gene_plot <- AnnotationPlot(
    object = atac, 
    region = peak_plot_expanded)
  peak_plot <- PeakPlot(
    object = atac,
    region = peak_plot_expanded, 
    peaks = Signac::StringToGRanges(peak_plot, sep = c("-", "-")))
  p <- CombineTracks(
    plotlist = list(cov_plot, peak_plot, gene_plot),
    heights = c(10, 1, 1),
  )
  ## write output
  file2write <- paste0(dir_out, gene_tmp, "_", chr, "_", st, "_", en, ".png")
  png(file2write, width = 1000, height = 1000, res = 150)
  print(p)
  dev.off()
}


