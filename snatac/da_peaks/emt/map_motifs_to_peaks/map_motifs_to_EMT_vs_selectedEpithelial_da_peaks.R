# Yige Wu @WashU Sep 2021
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
## set run id
version_tmp <- 3
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input daps
peaks_anno_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/emt/annotate_peaks/annotate_EMT_vs_selectedEpithelial_diff_peaks_to_genes/20210927.v1/ccRCC_vs_PT_DAPs.Annotated.20210927.v1.tsv")
## input peak-to-motif table
peak2motif_df <- fread(data.table = F, input = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/Signac.1.0.0/3.Merge_snATAC/Merge.SelectPeaks.v.20210706/peaks/Motifs_matched.28_snATAC_merged.object.20210827.tsv")

# overlap -----------------------------------------------------------------
peaks_anno_df <- merge(x = peaks_anno_df, y = peak2motif_df, by.x = c("peak"), by.y = c("Peak"), all.x = T)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "EMT_vs_selectedEpithelial_diff_peaks_to_motifs.", run_id, ".tsv")
write.table(x = peaks_anno_df, file = file2write, sep = "\t", row.names = F, quote = F)
