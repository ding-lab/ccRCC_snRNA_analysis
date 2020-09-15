# Yige Wu @WashU Sep 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
dir_input <- paste0("../ccRCC_snATAC/Resources/snATAC_Processed_Data/4.Cell_type_markers/4.DA_peaks/Major_cell_groups/out")

# get paths to the files of interest -------------------------------------------
paths_process <- list.files(path = dir_input, full.names = T)
paths_process

# input file and unite ----------------------------------------------------
peaks_df <- NULL
for (path_process in paths_process) {
  peaks_df_tmp <- fread(data.table = F, input = path_process)
  # peaks_df_tmp <- read.table(file = path_process)
  # peaks_df_tmp <- as.data.frame(peaks_df_tmp)
  rownames(peaks_df_tmp)
  ## get the cell types
  path_process
  filename_tmp <- str_split(string = path_process, pattern = "\\/")[[1]]
  filename_tmp <- filename_tmp[length(filename_tmp)]
  celltypename <- str_split(string = filename_tmp, pattern = "_")[[1]][1]
  celltypename
  ## add the cell type
  peaks_df_tmp <- peaks_df_tmp %>%
    # dplyr::select(p_val) %>%
    dplyr::mutate(Cell_type.filename = celltypename)
  ## unite
  # rownames_original <- c(rownames(peaks_df), rownames(peaks_df_tmp))
  # rownames_uniq <- make.names(names = rownames_original, unique = T)
  # rownames(peaks_df_tmp) <- rownames_uniq[(ifelse(is.null(peaks_df), 0, nrow(peaks_df))+1):length(rownames_original)]
  peaks_df <- rbind(peaks_df, peaks_df_tmp)
}
## add columns
peaks_df <- peaks_df %>%
  mutate(pct_diff = (pct.1 - pct.2))

# write output ------------------------------------------------------------
file_write <- paste0(dir_out, "DA_peaks.chromvar.MergedObj.byCell_group4.", run_id, ".tsv")
write.table(x = peaks_df, file = file_write, quote = F, sep = "\t", row.names = F)


