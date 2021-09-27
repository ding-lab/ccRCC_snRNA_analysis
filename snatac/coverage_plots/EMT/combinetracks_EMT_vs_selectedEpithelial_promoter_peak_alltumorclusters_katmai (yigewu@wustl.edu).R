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

# preprocess ATAC object --------------------------------------------------
atac@meta.data$cell_group <- mapvalues(x = rownames(atac), from = barcode2cluster_df$sample_barcode, to = as.vector(barcode2cluster_df$cell_group))
atac@meta.data$cell_group[atac@meta.data$cell_group == rownames(atac)] <- "other"
Idents(atac)=atac$cell_group
Idents(atac_subset)=factor(atac_subset$Piece_ID, levels=c(pieceids_selected,'C3L-00088-N','C3N-01200-N', "C3L-00079-N", "C3N-00242-N"))

# process coordinates ------------------------------------------------------------
chr=strsplit(x = peak_plot, split = "\\-")[[1]][1]
st=strsplit(x = peak_plot, split = "\\-")[[1]][2]; st = as.numeric(st)
en=strsplit(x = peak_plot, split = "\\-")[[1]][3]; en = as.numeric(en)
new_st=st-1000
new_en=en+1000
peak_plot_expanded=paste(chr,new_st,new_en,sep='-')
## process motif coordinates
motif_coord <- peak2motif_df$motif_coord[peak2motif_df$Peak == peak_plot & peak2motif_df$motif.name == motif_plot & peak2motif_df$Peak_Type == "Promoter"]; motif_coord <- unique(motif_coord)
## change atac ident
# print(head(atac@meta.data))

# make colors -------------------------------------------------------------
## make colors
color_tumorcell <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")[4]
color_pt <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")[1]
colors_celltype <- c(rep(x = color_tumorcell, 24), rep(x = color_pt, 6))
names(colors_celltype) <- c(peak2fcs_long_tmp_df$pieceid, sampleids_nat)

# plot --------------------------------------------------------------------
cov_plot= Signac::CoveragePlot(
  object = atac_subset,
  region = peak_plot_expanded,
  annotation = F, 
  peaks = F,
  links=FALSE)
cov_plot <- cov_plot + scale_fill_manual(values =  colors_celltype)
print("Finished cov_plot")

peak_plot_obj <- Signac::PeakPlot(
  object = atac_subset,
  region = peak_plot_expanded, 
  peaks = StringToGRanges(peak_plot, sep = c("-", "-")))
print("Finished peak plot")

gene_plot <- Signac::AnnotationPlot(
  object = atac_subset,
  region = peak_plot_expanded)

p <- Signac::CombineTracks(
  plotlist = list(cov_plot, peak_plot_obj, gene_plot),
  heights = c(7, 0.5, 0.5, 1))
print("Finished CombineTracks")

## write output
# file2write <- paste0(dir_out, gsub(x = peak_plot, pattern = "\\-", replacement = "_"), ".", motif_plot, ".png")
# png(file2write, width = 1000, height = 800, res = 150)
# print(p)
# dev.off()
file2write <- paste0(dir_out, gsub(x = peak_plot, pattern = "\\-", replacement = "_"), ".", motif_plot, ".pdf")
pdf(file2write, width = 6, height = 20)
print(p)
dev.off()
