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
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input -------------------------------------------------------------------
dirs_input <- c("/diskmnt/primary/ccRCC_snRNA/", "/diskmnt/primary/ccRCC_snATAC/")
fastq_summary_df <- NULL
fastq_summary_list <- list()
for (dir_input in dirs_input) {
  files_input <- list.files(path = dir_input, recursive = T)
  files_input <- files_input[grepl(pattern = "Samplemap\\.csv", x = files_input) & !grepl(pattern = "md5", x = files_input)]
  # print(files_input)
  datatype_tmp <- ifelse(dir_input == "/diskmnt/primary/ccRCC_snRNA/", "snRNA", "snATAC")
  for (file_tmp in files_input) {
    print(file_tmp)
    df_tmp <- fread(data.table = F, input = paste0(dir_input, file_tmp))
    df_tmp$data_type <- datatype_tmp
    df_tmp$file_dir <- str_split_fixed(string = file_tmp, pattern = "\\/", n = 2)[,1]
    print(head(df_tmp))
    colname_lane <- colnames(df_tmp)[grepl(pattern = "Lane", x = colnames(df_tmp))]
    colname_samplename <- colnames(df_tmp)[grepl(pattern = "Sample", x = colnames(df_tmp))]
    
    df_tmp2 <- df_tmp[c("Flow Cell ID", "Index Sequence", "Completion Date", colname_samplename, colname_lane, "Library Name", "data_type", "file_dir")]
    colnames(df_tmp2) <- c("Flow Cell ID", "Index Sequence", "Completion Date", "Sample Name", "Lane Number", "Library Name", "data_type", "file_dir")
    # df_tmp2 <- df_tmp[c("Flow Cell ID", "Index Sequence", "Library Name", "Completion Date", "data_type", "file_dir")]
    df_tmp2[, "Completion Date"] <- as.character(df_tmp2[, "Completion Date"])
    print(head(df_tmp2))
    
    fastq_summary_list[[file_tmp]] <- df_tmp
    fastq_summary_df <- rbind(fastq_summary_df, df_tmp2)
  }
}


# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "ccRCC.snRNA.snATAC.fastq.summary.", run_id, ".tsv")
write.table(x = fastq_summary_df, file = file2write, sep = "\t", row.names = F, quote = F)
file2write <- paste0(dir_out, "ccRCC.snRNA.snATAC.fastq.summary.", run_id, ".RSD")
saveRDS(object = fastq_summary_list, file = file2write, compress = T)