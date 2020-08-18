# Yige Wu @WashU June 2020

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

# input files ------------------------------------------------------
dir_files <- "./Resources/Analysis_Results/findmarkers/tumor_subclusters/"
filepaths <- list.files(dir_files, full.names = T, recursive = T)
filepaths
filepaths2process <- filepaths[grepl(x = filepaths, pattern = "findallmarkers_C3L-00010_correlated_tumorcells_newcluster0vs")]
filepaths2process
deg_all_df <- NULL
for (filepath in filepaths2process) {
  deg_tmp_df <- fread(data.table = F, input = filepath)
  deg_all_df <- rbind(deg_all_df, deg_tmp_df)
}

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "findallmarkers_C3L-00010_correlated_tumorcells_newcluster0vsothers.", run_id, ".tsv")
write.table(x = deg_all_df, file = file2write, quote = F, row.names = F, sep = "\t")
