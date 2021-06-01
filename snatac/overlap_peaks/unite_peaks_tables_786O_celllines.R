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
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input files -------------------------------------------------------------
samples <- c("BAP1-20210412-786-O", "control-20210412-786-O")
all_peaks=NULL
for (sample in samples){
  peaks=read.table(paste("/diskmnt/Projects/ccRCC_scratch/ccRCC_snATAC/Resources/snATAC_Processed_Data/Signac.1.0.0/1.Create_rds/out/",
                         sample,"/recentered_final.filtered",sample,".tsv",sep=""),
                   sep='\t',header=TRUE)
  peaks$Sample=sample
  all_peaks=rbind(all_peaks,peaks)
}

# write output ------------------------------------------------------------
write.table(x = all_peaks, file = paste0(dir_out, 'MACS2_peaks.BAP1_CellLines.BySample.',run_id, '.tsv'),
            sep='\t',quote=FALSE,row.names=FALSE)


