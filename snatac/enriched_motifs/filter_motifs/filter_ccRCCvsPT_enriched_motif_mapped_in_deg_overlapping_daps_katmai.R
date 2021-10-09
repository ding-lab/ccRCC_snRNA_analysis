# Yige Wu @WashU May 2021

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
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)


# input dependencies ------------------------------------------------------
## input motif mapped
motifs_all_df <- fread(data.table = F, input = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/Signac.1.0.0/3.Merge_snATAC/Merge.SelectPeaks.v.20210706/peaks/Motifs_matched.28_snATAC_merged.object.20210827.tsv")
# > head(motifs_all_df)
# group_name                     Peak strand    score              motif_coord
# 1   MA0003.4 chr5-135603208-135603708      + 11.45700 chr5-135603517-135603530
# 2   MA0003.4  chr19-47005145-47005645      - 16.12478  chr19-47005391-47005404
# 3   MA0003.4  chr10-13562850-13563350      + 12.12968  chr10-13563006-13563019
# 4   MA0003.4   chr7-47500038-47500538      + 11.63274   chr7-47500245-47500258
# 5   MA0003.4  chr14-91573795-91574295      + 13.75193  chr14-91574100-91574113
# 6   MA0003.4 chrX-130066889-130067389      + 13.85272 chrX-130067008-130067021
# motif.name
# 1     TFAP2A
# 2     TFAP2A
# 3     TFAP2A
# 4     TFAP2A
# 5     TFAP2A
# 6     TFAP2A
## input peak annotation
peaks_anno_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/annotate_peaks/annotate_ccrcc_vs_pt_promoter_daps/20211004.v1/ccRCC_vs_PTDAPs.Annotated.Promoter.20211004.v1.tsv")
## specify the motifs to filter
motifs_plot <- c("NFKB2", "NFKB1", "HIF1A", "ARNT::HIF1A", "RBPJ", "MXI1", "ZNF75D", "HSF2", "NEUROD1", "SREBF2", "NEUROG2(var.2)", "TBXT", "REL", "RELA", "KLF9")

# filter co-accessible peaks ----------------------------------------------
motifs_filtered_df <- motifs_all_df %>%
  filter(motif.name %in% motifs_plot)

dap2motif_df <- merge(x = peaks_anno_df, y = motifs_filtered_df, by.x = c("peak"), by.y = c("Peak"))

# write outupt ------------------------------------------------------------
file2write <- paste0(dir_out, "EnrichedMotifs2DAPs.tsv")
write.table(x = dap2motif_df, file = file2write, quote = F, sep = "\t", row.names = F)
