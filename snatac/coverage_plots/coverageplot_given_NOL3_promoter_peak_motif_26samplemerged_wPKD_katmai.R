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
atac=readRDS(paste('/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/Signac.1.0.0/3.Merge_snATAC/Merge.PKDSamples.v.20210609/26_ccRCC_4PKD_snATAC_ccRCC_peaks.20210611.rds.gz',sep=''))
Idents(atac)=atac$Piece_ID
## input motif-peak mapping result
peak2motif_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Motifs_Mapped_to_Peaks/Motifs_matched.DEG_associated_Peaks.Motif_annotation.20210517.v1.tsv")
## input peak fold changes
peak2fcs_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/ccRCC_Specific/DA_peaks_Tumor_vs_PT_affected_byCNV_removed.tsv")
## specify parameters to plot
peak_plot <- c("chr16-67170198-67170698")
motif_plot <- "HIF1A"
topn_plot <- 3

# preprocess samples to show ----------------------------------------------
peak2fcs_tmp_df <- peak2fcs_df %>%
  filter(peak == peak_plot)
peak2fcs_long_tmp_df <- melt(data = peak2fcs_tmp_df, measure.vars = colnames(peak2fcs_tmp_df)[grepl(pattern = "avg_lnFC", x = colnames(peak2fcs_tmp_df))])
peak2fcs_long_tmp_df <- peak2fcs_long_tmp_df %>%
  arrange(desc(value)) %>%
  mutate(pieceid = str_split_fixed(string = variable, pattern = "_", n = 2)[,1])
pieceids_selected <- head(x = peak2fcs_long_tmp_df$pieceid, topn_plot)

# preprocess ATAC object --------------------------------------------------
head(atac@meta.data)
atac_subset=subset(atac,(cell_type.v20210611 %in% c('Tumor') & Piece_ID %in% pieceids_selected) | cell_type.v20210611=='PT' & Piece_ID %in% c('C3L-00088-N','C3N-01200-N', "K1103044", "K1301462", "K1301463FB", "K1900070_1FB"))

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
Idents(atac_subset)=factor(atac_subset$Piece_ID, levels=c(pieceids_selected, 'C3L-00088-N','C3N-01200-N', "K1103044", "K1301462", "K1301463FB", "K1900070_1FB"))

# plot --------------------------------------------------------------------
p=Signac::CoveragePlot(
  object = atac_subset,
  region = peak_plot_expanded, 
  annotation = TRUE,
  peaks = TRUE,
  ranges=StringToGRanges(motif_coord, sep = c("-", "-")),
  ranges.title = paste(motif_plot,"\nmotif",sep=''),
  links=FALSE)

## write output
# file2write <- paste0(dir_out, gsub(x = peak_plot, pattern = "\\-", replacement = "_"), ".", motif_plot, ".png")
# png(file2write, width = 1000, height = 800, res = 150)
# print(p)
# dev.off()
file2write <- paste0(dir_out, gsub(x = peak_plot, pattern = "\\-", replacement = "_"), ".", motif_plot, ".pdf")
pdf(file2write, width = 6, height = 6)
print(p)
dev.off()

