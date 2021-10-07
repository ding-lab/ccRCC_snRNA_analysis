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
library(ggplot2)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the merged object
atac=readRDS(paste('/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/Signac.1.0.0/3.Merge_snATAC/Merge.SelectPeaks.v.20210526/',
                   '26_ccRCC_snATAC.selectedPeaks.chromvar.cicero.v4.20210607.rds',sep=''))
## input peak fold changes
peak2fcs_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Differential_Peaks/BAP1_Specific/DOWN_BAP1_specific_peaks.Filtered.CNV_corrected.Annotated.20210610.tsv")
## input mutation group
mut_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/mutation/annotate_cptac_sample_by_pbrm1_bap1_mutation/20210412.v1/PBRM1_BAP1_Mutation_Status_By_Case.20210412.v1.tsv")
## specify parameters to plot
peak_plot <- c("chr16-66955443-66955943")

# preprocess ATAC object --------------------------------------------------
pieceids_tumor_selected <- c("C3L-00908-T1", "C3L-00416-T2", ## BAP1&PBRM1 mutants
                             "C3L-01313-T1", "C3N-01200-T1", "C3N-00317-T1", "C3N-00437-T1", ## BAP1 mutants
                             "C3N-00733-T1", "C3L-00610-T1", "C3L-00079-T1", "C3N-00242-T1", "C3L-01302-T1", "C3N-01213-T1", "C3L-00004-T1", "C3L-00790-T1", "C3L-00583-T1",
                             "C3L-00917-T1", "C3L-00088-T1", "C3L-00088-T2", "C3L-00448-T1", "C3L-00096-T1", "C3L-00010-T1", "C3N-00495-T1", "C3L-00026-T1")
pieceids_nat_selected <- c('C3L-00088-N', "C3L-00079-N", "C3N-00242-N", "C3N-01200-N")

atac_subset=subset(atac,(cell_type %in% c('Tumor') & Piece_ID %in% pieceids_tumor_selected) | (cell_type=='PT' & Piece_ID %in% pieceids_nat_selected))
Idents(atac_subset)=factor(atac_subset$Piece_ID,levels=c(pieceids_nat_selected, pieceids_tumor_selected))

# process coordinates ------------------------------------------------------------
chr=strsplit(x = peak_plot, split = "\\-")[[1]][1]
st=strsplit(x = peak_plot, split = "\\-")[[1]][2]; st = as.numeric(st)
en=strsplit(x = peak_plot, split = "\\-")[[1]][3]; en = as.numeric(en)
new_st=st-1000
new_en=en+1000
peak_plot_expanded=paste(chr,new_st,new_en,sep='-')

# make colors ------------------------------------------------------------
## make colors
colors_tumorgroup_sim <- c(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[c(1, 2, 4)], RColorBrewer::brewer.pal(n = 8, name = "Set1")[c(4, 1)])
names(colors_tumorgroup_sim) <- c("PT",  "PBRM1 mutated", "Non-mutants", "BAP1 mutated", "Both mutated")
tumorgroup_vec <- mapvalues(x = str_split_fixed(string = pieceids_tumor_selected, pattern = "\\-T", n = 2)[,1], from = mut_df$Case, to = as.vector(mut_df$mutation_category_sim))
tumorgroup_vec <- c(rep("PT", length(pieceids_nat_selected)), tumorgroup_vec)
colors_tumorgroup <- colors_tumorgroup_sim[tumorgroup_vec]
names(colors_tumorgroup) <- c(pieceids_nat_selected, pieceids_tumor_selected)
print(colors_tumorgroup)

# plot --------------------------------------------------------------------
p=CoveragePlot(
  object = atac_subset,
  region = peak_plot_expanded, 
  annotation = TRUE,
  peaks = TRUE,
  links=FALSE, heights = c(10, 1, 0.5))
p <- p + scale_fill_manual(values =  colors_celltype)

## write output
file2write <- paste0(dir_out, gsub(x = peak_plot, pattern = "\\-", replacement = "_"), ".pdf")
pdf(file2write, width = 6, height = 10, useDingbats = F)
print(p)
dev.off()

