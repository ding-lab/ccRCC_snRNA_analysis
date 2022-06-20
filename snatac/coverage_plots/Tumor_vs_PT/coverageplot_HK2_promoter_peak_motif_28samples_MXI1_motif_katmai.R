# Yige Wu @WashU May 2021
## source activate signac

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
setwd(dir_base)
## load libraries
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "Signac",
  "Seurat",
  "ggplot2"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    print(paste0(pkg_name_tmp, "is being installed!"))
    BiocManager::install(pkgs = pkg_name_tmp, update = F)
    install.packages(pkg_name_tmp, dependencies = T)
  }
  print(paste0(pkg_name_tmp, " is installed!"))
  library(package = pkg_name_tmp, character.only = T)
}
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the merged object
atac=readRDS(paste('/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/Signac.1.0.0/3.Merge_snATAC/Merge.SelectPeaks.v.20210706/28_ccRCC_snATAC.selectedPeaks.chromvar.cicero.v3.20210725.rds',sep=''))
Idents(atac)=atac$Piece_ID
## input motif-peak mapping result
# peak2motif_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Motifs_Mapped_to_Peaks/Motifs_matched.DEG_associated_Peaks.Motif_annotation.20210517.v1.tsv")
peak2motif_df <- fread(data.table = F, input = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/Signac.1.0.0/3.Merge_snATAC/Merge.SelectPeaks.v.20210706/peaks/Motifs_matched.28_snATAC_merged.object.20210827.tsv")
## input peak fold changes
peak2fcs_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/ccRCC_Specific/DA_peaks_Tumor_vs_PT_affected_byCNV_removed.tsv")

# specify parameters ------------------------------------------------------
motifs_plot <- c("MXI1")
motifs_plot <- c("MA1108.2")
peak_plot <- c("chr2-74876177-74876677")
topn_plot <- 24
motif_coord <- peak2motif_df$motif_coord[(peak2motif_df$Peak %in% peak_plot) & (peak2motif_df$group_name %in% motifs_plot)]
motif_coord <- unique(motif_coord)
motif_coord <- sort(motif_coord); motif_coord 
# [1] "chr2-74834043-74834058" "chr2-74856215-74856224" "chr2-74876328-74876337"
## KLF9, MXI1, MXI1


# preprocess samples to show ----------------------------------------------
pieceids_tumor_selected <- c("C3L-00448-T1", "C3L-01302-T1","C3L-00088-T1", "C3N-00242-T1", "C3L-00790-T1", "C3L-00088-T2", 
                             "C3L-01313-T1", "C3L-00917-T1", "C3N-01200-T1", "C3L-00610-T1", "C3N-01213-T1", "C3L-00079-T1", 
                             "C3L-00583-T1",  "C3N-00733-T1", "C3N-00317-T1", "C3L-00908-T1", "C3L-00026-T1", "C3L-01287-T1", 
                             "C3L-00004-T1", "C3L-00416-T2", "C3L-00096-T1", "C3N-00437-T1", "C3L-00010-T1", "C3N-00495-T1")
pieceids_nat_selected <- c("C3N-00242-N", 'C3L-00088-N', "C3L-00079-N", 'C3N-01200-N')

# preprocess ATAC object --------------------------------------------------
head(atac@meta.data)
atac_subset=subset(atac,(cell_type %in% c('Tumor') & (Piece_ID %in% pieceids_tumor_selected)) | cell_type=='PT' & Piece_ID %in% pieceids_nat_selected)
## make colors

color_tumorcell <- RColorBrewer::brewer.pal(n = 9, name = "Dark2")[4]
color_pt <- RColorBrewer::brewer.pal(n = 9, name = "Dark2")[1]
colors_celltype <- c(rep(x = color_tumorcell, length(pieceids_tumor_selected)), rep(x = color_pt, length(pieceids_nat_selected)))
names(colors_celltype) <- c(pieceids_tumor_selected, pieceids_nat_selected)

# process coordinates ------------------------------------------------------------
chr=strsplit(x = peak_plot, split = "\\-")[[1]][1]
st=strsplit(x = peak_plot, split = "\\-")[[1]][2]; st = as.numeric(st)
en=strsplit(x = peak_plot, split = "\\-")[[1]][3]; en = as.numeric(en)
new_st=st-1000
new_en=en+1000
peak_plot_expanded=paste(chr,new_st,new_en,sep='-')
# peak_plot_expanded <- "chr2-74833572-74876777"
## change atac ident
# print(head(atac@meta.data))
Idents(atac_subset)=factor(atac_subset$Piece_ID, levels=c(pieceids_tumor_selected, pieceids_nat_selected))

# plot --------------------------------------------------------------------
cov_plot= Signac::CoveragePlot(
  object = atac_subset,
  region = peak_plot_expanded,
  annotation = F, 
  peaks = F,
  links=FALSE)
cov_plot <- cov_plot + scale_fill_manual(values =  colors_celltype)
print("Finished cov_plot")

peakplot_obj <- Signac::PeakPlot(
  object = atac_subset,
  region = peak_plot_expanded, 
  peaks = StringToGRanges(peak_plot, sep = c("-", "-")))
print("Finished peak plot")

motifplot_obj <- Signac::PeakPlot(
  object = atac_subset,
  region = peak_plot_expanded, 
  peaks = StringToGRanges(motif_coord, sep = c("-", "-")))
print("Finished motif plot")

gene_plot <- Signac::AnnotationPlot(
  object = atac_subset,
  region = peak_plot_expanded)

p <- Signac::CombineTracks(
  plotlist = list(cov_plot, peakplot_obj, motifplot_obj, gene_plot),
  heights = c(8, 0.5, 0.2, 1))
print("Finished CombineTracks")

## write output
# file2write <- paste0(dir_out, gsub(x = peak_plot, pattern = "\\-", replacement = "_"), ".", motif_plot, ".png")
# png(file2write, width = 1000, height = 800, res = 150)
# print(p)
# dev.off()
file2write <- paste0(dir_out, gsub(x = peak_plot[1], pattern = "\\-", replacement = "_"), ".", paste0(motifs_plot, collapse = "_"), ".pdf")
pdf(file2write, width = 6, height = 10, useDingbats = F)
print(p)
dev.off()

