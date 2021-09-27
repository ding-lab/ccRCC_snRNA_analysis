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
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
library(Signac)
library(ggplot2)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the merged object
atac=readRDS(paste('/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/Signac.1.0.0/3.Merge_snATAC/Merge.SelectPeaks.v.20210706/28_ccRCC_snATAC.selectedPeaks.chromvar.cicero.v3.20210725.rds',sep=''))
## input barcode to tumor clusters
barcode2cluster_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Barcode_Annotation/RNA_clusters_Tumor_PT_clusters.snATAC.20210917.tsv")
# barcode2cluster_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Barcode_Annotation/28_ccRCC_snATAC_ManualReviwed.v2.20210709.tsv")
## input the peaks to plot
peaks2degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/emt/overlap_degs/overlap_EMT_vs_selectedEpithelial_diff_promoter_peaks_with_degs/20210927.v1/EMT_vs_selectedEpithelialClusters.20210927.v1.tsv")

# preprocess ATAC object --------------------------------------------------
atac@meta.data$cell_group <- mapvalues(x = rownames(atac), from = barcode2cluster_df$sample_barcode, to = as.vector(barcode2cluster_df$cell_group))
atac@meta.data$cell_group[atac@meta.data$cell_group == rownames(atac)] <- "other"
Idents(atac)=atac$cell_group
atac_subset=subset(atac,(cell_group %in% c("C3L-00079-T1_C4", "C3L-01302-T1_C1",
                                           "C3L-00088-T2_C1", "C3N-00733-T1_C1", "C3L-00416-T2_C1", "C3L-00010-T1_C1", "C3L-00088-T1_C1")))

# process peaks -----------------------------------------------------------
plotdata_df <- peaks2degs_df %>%
  filter(!is.na(avg_log2FC.snATAC) & !is.na(avg_log2FC.snRNA)) %>%
  mutate(highlight = (avg_log2FC.snATAC >= 1 & avg_log2FC.snRNA >= 1) | (avg_log2FC.snATAC <= -1 & avg_log2FC.snRNA <= -1) | (Gene %in% c("VIM", "FN1", "CDH2", "WNT5B"))) %>%
  filter(highlight == T) %>%
  select(avg_log2FC.snATAC, avg_log2FC.snRNA, Gene, peak) %>%
  unique()

for (peak_plot in unique(plotdata_df$peak)) {
  # process coordinates ------------------------------------------------------------
  chr=strsplit(x = peak_plot, split = "\\-")[[1]][1]
  st=strsplit(x = peak_plot, split = "\\-")[[1]][2]; st = as.numeric(st)
  en=strsplit(x = peak_plot, split = "\\-")[[1]][3]; en = as.numeric(en)
  new_st=st-1000
  new_en=en+1000
  peak_plot_expanded=paste(chr,new_st,new_en,sep='-')
  gene_plot <- plotdata_df$Gene[plotdata_df$peak == peak_plot_expanded]
  
  # plot --------------------------------------------------------------------
  cov_obj= Signac::CoveragePlot(
    object = atac_subset,
    region = peak_plot_expanded,
    annotation = F, 
    peaks = F,
    links=FALSE)
  # cov_obj <- cov_obj + scale_fill_manual(values =  colors_celltype)
  print("Finished cov_plot")
  
  peak_plot_obj <- Signac::PeakPlot(
    object = atac_subset,
    region = peak_plot_expanded, 
    peaks = StringToGRanges(peak_plot, sep = c("-", "-")))
  print("Finished peak plot")
  
  gene_plot_obj <- Signac::AnnotationPlot(
    object = atac_subset,
    region = peak_plot_expanded)
  
  p <- Signac::CombineTracks(
    plotlist = list(cov_obj, peak_plot_obj, gene_plot_obj),
    heights = c(7, 0.5, 0.5, 1))
  print("Finished CombineTracks")
  
  ## write output
  # file2write <- paste0(dir_out, gsub(x = peak_plot, pattern = "\\-", replacement = "_"), ".", motif_plot, ".png")
  # png(file2write, width = 1000, height = 800, res = 150)
  # print(p)
  # dev.off()
  file2write <- paste0(dir_out, gene_plot, "_", gsub(x = peak_plot, pattern = "\\-", replacement = "_"), ".pdf")
  pdf(file2write, width = 6, height = 8)
  print(p)
  dev.off()
}
