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
library(Signac)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 2
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
peak2fcs_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/BAP1_Specific/DOWN_BAP1_specific_peaks.Filtered.CNV_corrected.Annotated.20210610.tsv")
# peaks_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/bap1/overlap_bap1_specific_enhancer_promoter_peaks_with_degs/20210615.v1/BAP1_DAP2DEG.20210615.v1.tsv")
## specify parameters to plot
peak_plot <- c("chr16-66955443-66961444")
topn_plot <- 3

# preprocess samples to show ----------------------------------------------
peak2fcs_tmp_df <- peak2fcs_df %>%
  filter(peak == peak_plot)
peak2fcs_long_tmp_df <- melt(data = peak2fcs_tmp_df, measure.vars = colnames(peak2fcs_tmp_df)[grepl(pattern = "_Signif_avg_lnFC", x = colnames(peak2fcs_tmp_df))])
peak2fcs_long_tmp_df <- peak2fcs_long_tmp_df %>%
  arrange(value) %>%
  mutate(pieceid = str_split_fixed(string = variable, pattern = "_", n = 2)[,1])
pieceids_selected <- head(x = peak2fcs_long_tmp_df$pieceid, topn_plot)
pieceids_selected <- c(pieceids_selected, "C3N-00733-T1", "C3N-00242-T1", "C3L-00917-T1", "C3L-00026-T1")

# preprocess ATAC object --------------------------------------------------
atac_subset=subset(atac,(cell_type %in% c('Tumor') & Piece_ID %in% pieceids_selected) | cell_type=='PT' & Piece_ID %in% c('C3L-00088-N','C3N-01200-N'))

# process coordinates ------------------------------------------------------------
chr=strsplit(x = peak_plot, split = "\\-")[[1]][1]
st=strsplit(x = peak_plot, split = "\\-")[[1]][2]; st = as.numeric(st)
en=strsplit(x = peak_plot, split = "\\-")[[1]][3]; en = as.numeric(en)
new_st=st-500
new_en=en+500
peak_plot_expanded=paste(chr,new_st,new_en,sep='-')
## change atac ident
# print(head(atac@meta.data))
Idents(atac_subset)=factor(atac_subset$Piece_ID,levels=c(pieceids_selected, 'C3L-00088-N','C3N-01200-N'))

# plot --------------------------------------------------------------------
p=Signac::CoveragePlot(
  object = atac_subset,
  region = peak_plot_expanded, 
  annotation = TRUE,
  peaks = TRUE,
  links=FALSE)

## write output
file2write <- paste0(dir_out, gsub(x = peak_plot, pattern = "\\-", replacement = "_"), ".pdf")
pdf(file2write, width = 6, height = 7, useDingbats = F)
print(p)
dev.off()

