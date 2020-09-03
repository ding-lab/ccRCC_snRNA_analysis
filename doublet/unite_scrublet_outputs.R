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
## input meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")
## get file paths
dir_scrublet_out <- "./Resources/Analysis_Results/scrublet/run20200902_adj_cutoff/"
paths_file <- list.files(dir_scrublet_out, full.names = T)
paths_file
paths_file_process <- paths_file[grepl(pattern = "output_table.csv", x = paths_file)]
paths_file_process

# unite -------------------------------------------------------------------
scrublet_united_df <- NULL
for (path_file in paths_file_process) {
  scrublet_df <- fread(data.table = F, input = path_file)
  ## get aliquot id
  id_aliquot <- str_split_fixed(string = path_file, pattern = dir_scrublet_out, n = 2)[,2]
  id_aliquot
  id_aliquot <- str_split_fixed(string = id_aliquot, pattern = "_|\\/", n = 3)[,2]
  id_aliquot_wu <- idmetadata_df$Aliquot.snRNA.WU[idmetadata_df$Aliquot.snRNA == id_aliquot]
  ## add aliquot id and reable aliquot id
  scrublet_df$Aliquot <- id_aliquot
  scrublet_df$Aliquot_WU <- id_aliquot_wu
  ## unite
  scrublet_united_df <- rbind(scrublet_united_df, scrublet_df)
}
# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "scrublet.", "run20200902_adj_cutoff.", "united_outputs", ".tsv")
write.table(x = scrublet_united_df, file = file2write, quote = F, sep = "\t", row.names = F)
