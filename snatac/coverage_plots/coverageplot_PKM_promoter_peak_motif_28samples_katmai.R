# Yige Wu @WashU May 2021
## source activate ccrcc_snrna
## 2021-08-27 right now the RBPJ motif is not there

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
Idents(atac)=atac$Piece_ID
## input cell type annotation
# barcode2celltype_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Barcode_Annotation/28_ccRCC_snATAC_ManualReviwed.v2.20210709.tsv")
## input motif-peak mapping result
peak2motif_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Motifs_Mapped_to_Peaks/Motifs_matched.DEG_associated_Peaks.Motif_annotation.20210517.v1.tsv")
## input peak fold changes
# peak2fcs_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/ccRCC_Specific/DA_peaks_Tumor_vs_PT_affected_byCNV_removed.tsv")
peak2fcs_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/ccRCC_Specific/UP_Tumor_vsPT.Filtered.CNV_corrected.Annotated.20210811.tsv")
## specify parameters to plot
peak_plot <- c("chr15-72228266-72228766")
motif_plot <- "RBPJ"
topn_plot <- 24

# preprocess samples to show ----------------------------------------------
peak2fcs_tmp_df <- peak2fcs_df %>%
  filter(peak == peak_plot)
peak2fcs_long_tmp_df <- melt(data = peak2fcs_tmp_df, measure.vars = colnames(peak2fcs_tmp_df)[grepl(pattern = "_Signif_avg_lnFC", x = colnames(peak2fcs_tmp_df))])
peak2fcs_long_tmp_df <- peak2fcs_long_tmp_df %>%
  arrange(desc(value)) %>%
  mutate(pieceid = str_split_fixed(string = variable, pattern = "_", n = 2)[,1])
pieceids_selected <- head(x = peak2fcs_long_tmp_df$pieceid, topn_plot)

# preprocess ATAC object --------------------------------------------------
head(atac@meta.data)
atac_subset=subset(atac,(cell_type %in% c('Tumor') & Piece_ID %in% pieceids_selected) | cell_type=='PT' & Piece_ID %in% c('C3L-00088-N','C3N-01200-N', "C3L-00079-N", "C3N-00242-N"))
## change atac ident
# print(head(atac@meta.data))
Idents(atac_subset)=factor(atac_subset$Piece_ID, levels=c(pieceids_selected,'C3L-00088-N','C3N-01200-N', "C3L-00079-N", "C3N-00242-N"))

# process coordinates ------------------------------------------------------------
chr=strsplit(x = peak_plot, split = "\\-")[[1]][1]
st=strsplit(x = peak_plot, split = "\\-")[[1]][2]; st = as.numeric(st)
en=strsplit(x = peak_plot, split = "\\-")[[1]][3]; en = as.numeric(en)
new_st=st-1000
new_en=en+1000
peak_plot_expanded=paste(chr,new_st,new_en,sep='-')
## process motif coordinates
motif_coord <- peak2motif_df$motif_coord[peak2motif_df$Peak == "chr15-72230151-72230651" & peak2motif_df$motif.name == motif_plot & peak2motif_df$Peak_Type == "Promoter"]; motif_coord <- unique(motif_coord)
## make colors
color_tumorcell <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")[4]
color_pt <- RColorBrewer::brewer.pal(n = 8, name        = "Dark2")[1]
colors_celltype <- c(rep(x = color_tumorcell, 24), rep(x = color_pt, 4))
names(colors_celltype) <- c(pieceids_selected, 'C3L-00088-N','C3N-01200-N', "C3L-00079-N", "C3N-00242-N")

# plot --------------------------------------------------------------------
p=Signac::CoveragePlot(
  object = atac_subset,
  region = peak_plot_expanded, 
  annotation = TRUE,
  peaks = TRUE,
  ranges=StringToGRanges(motif_coord, sep = c("-", "-")),
  ranges.title = paste(motif_plot,"\nmotif",sep=''),
  links=FALSE)
# p <- p + scale_fill_manual(values =  colors_celltype)

## write output
# file2write <- paste0(dir_out, gsub(x = peak_plot, pattern = "\\-", replacement = "_"), ".", motif_plot, ".png")
# png(file2write, width = 1000, height = 800, res = 150)
# print(p)
# dev.off()
file2write <- paste0(dir_out, gsub(x = peak_plot, pattern = "\\-", replacement = "_"), ".", motif_plot, ".pdf")
pdf(file2write, width = 6, height = 20)
print(p)
dev.off()

