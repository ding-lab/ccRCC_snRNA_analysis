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
# library(Signac)
library(data.table)
library(stringr)
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 6
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the merged object
# atac <- readRDS(file = "./Resources/Analysis_Results/snatac/merge_objects/merge_objects_786O_celllines/20210527.v1/786O_CellLines.Merged.20210527.v1.RDS")
atac <- readRDS(file = "./Resources/Analysis_Results/snatac/merge_objects/add_gene_activity_to_786O_celllines_merged_katmai/20210528.v1/786O_CellLines.Merged.20210528.v1.RDS")
print("Finished readRDS")
## specify the peak to plot
peaks_celline_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/bap1/overlap_bap1_associated_peaks_with_cellline_peaks/20210611.v1/Consistent.Promoter.BAP1_associated_Peaks.HumanTumorTissue_CellLine.20210611.v1.tsv")

for (i in 1:nrow(peaks_celline_df)) {
# for (i in 1) {
  peak_cellline_tmp <- peaks_celline_df$peak.snATAC.cellline[i]
  peak_tumortissue_tmp <- peaks_celline_df$peak.snATAC.tumortissue[i]
  ## process peak range for cell line
  chr =str_split_fixed(string = peak_cellline_tmp, pattern = "\\-", n = 3)[,1]
  st_cl=str_split_fixed(string = peak_cellline_tmp, pattern = "\\-", n = 3)[,2]; st_cl = as.numeric(st_cl)
  en_cl=str_split_fixed(string = peak_cellline_tmp, pattern = "\\-", n = 3)[,3]; en_cl = as.numeric(en_cl)

  ## process peak range for cell line for tumor tissue
  st_tt=str_split_fixed(string = peak_tumortissue_tmp, pattern = "\\-", n = 3)[,2]; st_tt = as.numeric(st_tt)
  en_tt=str_split_fixed(string = peak_tumortissue_tmp, pattern = "\\-", n = 3)[,3]; en_tt = as.numeric(en_tt)
  
  new_st=(min(c(st_cl, st_tt)) - 1000)
  new_en=(max(c(en_cl, en_tt)) + 1000)
  peak_plot_expanded=paste(chr,new_st,new_en,sep='-')
  
  # plot --------------------------------------------------------------------
  cov_plot= Signac::CoveragePlot(
    object = atac, group.by = "Piece_ID", 
    region = peak_plot_expanded,
    annotation = F, 
    peaks = F,
    links=FALSE)
  print("Finished cov_plot")
  
  peak_cl_plot <- Signac::PeakPlot(
    object = atac,
    region = peak_plot_expanded, 
    peaks = Signac::StringToGRanges(peak_cellline_tmp, sep = c("-", "-")))
  print("Finished peak_plot for cell line")
  
  peak_tt_plot <- Signac::PeakPlot(
    object = atac,
    region = peak_plot_expanded, 
    peaks = Signac::StringToGRanges(peak_tumortissue_tmp, sep = c("-", "-")))
  print("Finished peak_plot for tumor tissue")
  
  gene_plot <- Signac::AnnotationPlot(
    object = atac,
    region = peak_plot_expanded)
  
  p <- Signac::CombineTracks(
    plotlist = list(cov_plot, peak_cl_plot, peak_tt_plot, gene_plot),
    heights = c(10, 1, 1, 1))
  print("Finished CombineTracks")
  
  ## write output
  gene_tmp <- peaks_celline_df$Gene[i]
  file2write <- paste0(dir_out, gene_tmp, 
                       ".CellLine_", chr, "_", st_cl, "_", en_cl,
                       ".Tumor_", chr, "_", st_tt, "_", en_tt,
                       ".png")
  png(file2write, width = 1000, height = 1000, res = 150)
  print(p)
  dev.off()
}


