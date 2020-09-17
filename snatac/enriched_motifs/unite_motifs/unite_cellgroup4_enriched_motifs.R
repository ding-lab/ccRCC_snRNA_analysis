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
dir_input <- paste0("../ccRCC_snATAC/Resources/snATAC_Processed_Data/4.Cell_type_markers/3.EnrichedMotifs/Major_cell_group/out")

# get paths to the files of interest -------------------------------------------
paths_process <- list.files(path = dir_input, full.names = T)
paths_process

# input file and unite ----------------------------------------------------
motifs_df <- NULL
for (path_process in paths_process) {
  motifs_df_tmp <- fread(data.table = F, input = path_process)
  # motifs_df_tmp <- read.table(file = path_process)
  # motifs_df_tmp <- as.data.frame(motifs_df_tmp)
  rownames(motifs_df_tmp)
  ## get the cell types
  path_process
  filename_tmp <- str_split(string = path_process, pattern = "\\/")[[1]]
  filename_tmp <- filename_tmp[length(filename_tmp)]
  celltypename <- str_split(string = filename_tmp, pattern = "_")[[1]][1]
  celltypename
  ## add the cell type
  motifs_df_tmp <- motifs_df_tmp %>%
    # dplyr::select(p_val) %>%
    dplyr::mutate(Cell_type.filename = celltypename)
  ## unite
  rownames_original <- c(rownames(motifs_df), rownames(motifs_df_tmp))
  rownames_uniq <- make.names(names = rownames_original, unique = T)
  rownames(motifs_df_tmp) <- rownames_uniq[(ifelse(is.null(motifs_df), 0, nrow(motifs_df))+1):length(rownames_original)]
  motifs_df <- rbind(motifs_df, motifs_df_tmp)
}
## add columns
motifs_df <- motifs_df %>%
  mutate(pct_diff = (pct.1 - pct.2))

# write output ------------------------------------------------------------
file_write <- paste0(dir_out, "Enriched_Motifs.chromvar.MergedObj.byCell_group4.", run_id, ".tsv")
write.table(x = motifs_df, file = file_write, quote = F, sep = "\t", row.names = F)


