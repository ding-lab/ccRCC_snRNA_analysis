# Yige Wu @WashU May 2021
## source activate ccrcc_snrna
## reference: https://satijalab.org/signac/articles/merging.html

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
library(doParallel)
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input united peaks
p1 <- read.table(file = "./Resources/Analysis_Results/snatac/overlap_peaks/unite_peaks_tables_786O_celllines/20210601.v1/MACS2_peaks.BAP1_CellLines.BySample.20210601.v1.tsv", sep='\t',header=TRUE)

# overlap -----------------------------------------------------------------
recentered_p=Signac::StringToGRanges(p1$new_peak, sep = c("-", "-"))
print(paste0("Finished StringToGRanges!"))

olap=as.data.frame(Signac::findOverlaps(recentered_p,recentered_p))
print(paste0("Finished findOverlaps!"))

olap1=olap[olap$queryHits!=olap$subjectHits,]

recentered_non_olap=p1[-olap1$queryHits,]

pairs=cbind(p1[olap1$queryHits,c(1:3,7)],olap1$queryHits,p1[olap1$subjectHits,c(1:3,7)],olap1$subjectHits)
colnames(pairs)=c('chr_1','st_1','en_1','score_1','row_1','chr_2','st_2','en_2','score_2','row_2')

pairs=pairs[pairs$score_1>=pairs$score_2,]
pairs=pairs[order(-pairs$score_1),]
pairs_all=pairs

# parallel ----------------------------------------------------------------
registerDoParallel(cores=28)

all_st=NULL
all_st<-foreach(chr_n=c(1:22,"X","Y")) %dopar% {
  chr=paste("chr",chr_n,sep='')
  pairs=pairs_all[pairs_all$chr_1==chr,]
  pairs=pairs[,c(4,5,9,10)]
  all_st_chr=NULL
  for (i in 1:nrow(pairs)){
    if (nrow(pairs)>0){
      p_del=pairs[pairs$row_1==pairs[1,2],]
      all_st_chr=rbind(all_st_chr,p1[rownames(p1)==pairs[1,2],])
      pairs=pairs[!(pairs$row_1 %in% c(p_del$row_1[1],p_del$row_2)),]
    }
    #    print(paste(chr,nrow(pairs),sep=' '))
  }
  return(all_st_chr)
}

all_st_f=NULL
for (i in 1:24){
  all_st_1=as.data.frame(all_st[[i]])
  all_st_1=all_st_1[!duplicated(all_st_1),]
  all_st_f=rbind(all_st_f,all_st_1)
  print(paste0("rbind ", i))
}
recentered_final=rbind(recentered_non_olap,all_st_f)

# write output ------------------------------------------------------------
write.table(recentered_final,paste(dir_out, 'recentered_final.filtered.', run_id,'.tsv',sep=''),sep='\t',
            quote=FALSE,row.names=FALSE)
write.table(recentered_non_olap,paste(dir_out, 'recentered_nonOverlapping.filtered.', run_id,'.tsv',sep=''),sep='\t',
            quote=FALSE,row.names=FALSE)
write.table(all_st_f,paste(dir_out, 'recentered_Overlapping.filtered.', run_id,'.tsv',sep=''),sep='\t',
            quote=FALSE,row.names=FALSE)