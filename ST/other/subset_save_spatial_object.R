## save object for HT283 and 293
## R version
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
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "Seurat",
  "patchwork",
  "tidyverse"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}

####################################### 
## Load and process
source('/diskmnt/Datasets/Spatial_Transcriptomics/Analysis/RCC/script/ccRCC_loadST_08092021_wo_RCTD.R')

## Subset HT282 and HT293
ht_ids = st_merge@meta.data %>% filter(str_detect(orig.ident, 'HT282|HT293')) %>% rownames
st_htan = subset(st_merge, cells = ht_ids)

# save --------------------------------------------------------------------
filename_write <- paste0(dir_out, "ccRCC_ST.HT282_HT293.", run_id, ".RDS")
saveRDS(object = st_htan, file = filename_write, compress = T)

