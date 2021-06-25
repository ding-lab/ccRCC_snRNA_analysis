# Yige Wu @WashU Jun 2021
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
version_tmp <- "1"
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the merged object
# atac=readRDS(paste('/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/Signac.1.0.0/3.Merge_snATAC/Merge.SelectPeaks.v.20210526/',
#                    '26_ccRCC_snATAC.selectedPeaks.chromvar.v3.20210602.rds',sep=''))
atac=readRDS(paste('/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/Signac.1.0.0/3.Merge_snATAC/Merge.SelectPeaks.v.20210503/',
                   '26_ccRCC_snATAC.selectedPeaks.chromvar.CICERo.v6.20210512.rds',sep=''))
Idents(atac)=atac$Piece_ID
## input peak fold changes
peaks_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/pbrm1/overlap_degs/overlap_pbrm1_vs_nonmutant_enhancer_promoter_peaks_with_degs/20210625.v1/PBRM1_vs_NonMutant_DAP2DEG.20210625.v1.tsv")

# preprocess ATAC object --------------------------------------------------
pieceids_selected <- c("C3N-00733-T1", "C3L-00079-T1", "C3L-01302-T1", "C3N-01213-T1", "C3L-0004-T1", "C3L-00790-T1", "C3L-00583-T1", "C3L-00610-T1", "C3N-00242-T1", ## PBRM1 mutants
                       "C3L-00908-T1", "C3L-00416-T2", ## PBRM1 & BAP1 mutants
                       "C3L-01313-T1", "C3N-01200-T1", "C3N-00317-T1", "C3N-00437-T1", ## BAP1 mutants
                       "C3L-00096-T1", "C3N-00495-T1", "C3L-00917-T1", "C3L-00088-T1", "C3L-00088-T2", "C3L-00448-T1", "C3L-00010-T1", "C3L-00026-T1")
atac_subset=subset(atac,(cell_type %in% c('Tumor') & Piece_ID %in% pieceids_selected) | cell_type=='PT' & Piece_ID %in% c('C3L-00088-N','C3N-01200-N'))
Idents(atac_subset)=factor(atac_subset$Piece_ID,levels=c(pieceids_selected, 'C3L-00088-N','C3N-01200-N'))

# make colors -------------------------------------------------------------
## make colors
color_tumorcell_pbrm1 <- RColorBrewer::brewer.pal(n = 9, name = "Set1")[5]
color_tumorcell_bap1 <- RColorBrewer::brewer.pal(n = 9, name = "Set1")[4]
color_tumorcell_nonmutnat <- RColorBrewer::brewer.pal(n = 9, name = "Dark2")[4]
color_pt <- RColorBrewer::brewer.pal(n = 9, name = "Dark2")[1]
colors_celltype <- c(rep(x = color_tumorcell_pbrm1, 9), 
                     rep(x = color_tumorcell_bap1, 6),
                     rep(x = color_tumorcell_nonmutnat, 8),
                     rep(x = color_pt, 2))
names(colors_celltype) <- c(pieceids_selected, 'C3L-00088-N','C3N-01200-N')

# process coordinates ------------------------------------------------------------
## specify parameters to plot
peaks_plot <- peaks_df$peak[peaks_df$peak2gene_type == "Promoter"]
peaks_plot <- peaks_df$peak

# for (peak_plot in c("chr16-66955443-66955943")) {
for (peak_plot in unique(peaks_plot)) {

  chr=strsplit(x = peak_plot, split = "\\-")[[1]][1]
  st=strsplit(x = peak_plot, split = "\\-")[[1]][2]; st = as.numeric(st)
  en=strsplit(x = peak_plot, split = "\\-")[[1]][3]; en = as.numeric(en)
  new_st=st-1000
  new_en=en+1000
  peak_plot_expanded=paste(chr,new_st,new_en,sep='-')
  gene_plot <- peaks_df$Gene[peaks_df$peak == peak_plot]
  peak2gene_type_plot <- peaks_df$peak2gene_type[peaks_df$peak == peak_plot]
  
  # plot --------------------------------------------------------------------
  cov_plot= Signac::CoveragePlot(
    object = atac_subset,
    region = peak_plot_expanded,
    annotation = F, 
    peaks = F,
    links=FALSE)
  cov_plot <- cov_plot + scale_fill_manual(values =  colors_celltype)
  print("Finished cov_plot")
  
  peak_obj <- Signac::PeakPlot(
    object = atac_subset,
    region = peak_plot_expanded, 
    peaks = StringToGRanges(peak_plot, sep = c("-", "-")))
  print("Finished peak plot")
  
  gene_obj <- Signac::AnnotationPlot(
    object = atac_subset,
    region = peak_plot_expanded)
  
  p <- Signac::CombineTracks(
    plotlist = list(cov_plot, peak_obj, gene_obj),
    heights = c(23, 0.5, 1.5))
  print("Finished CombineTracks")
  
  ## write output
  file2write <- paste0(dir_out, paste0(gene_plot, collapse = "_"), 
                       paste0(peak2gene_type_plot, collapse = "_"), 
                       gsub(x = peak_plot, pattern = "\\-", replacement = "_"), ".pdf")
  pdf(file2write, width = 6, height = 9, useDingbats = F)
  print(p)
  dev.off()
}



