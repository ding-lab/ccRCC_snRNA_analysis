# Yige Wu @WashU Apr 2021

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
# dir_deg_file_level1 <- "./Resources/snRNA_Processed_Data/Differentially_Expressed_Genes/Tumorcells_vs_PTcells/"
dir_deg_file_level1 <- "./Resources/snRNA_Processed_Data/Differentially_Expressed_Genes/Tumorcells_vs_PTcells/20210429.v2/"

# input files -------------------------------------------------------------
paths_deg_file <- list.files(path = dir_deg_file_level1, full.names = T)
paths_deg_file <- paths_deg_file[grepl(pattern = "tsv", x = paths_deg_file)]
deg_sup_df <- NULL
for (path_deg_file_tmp in paths_deg_file) {
  deg_df_tmp <- fread(data.table = F, input = path_deg_file_tmp)
  deg_sup_df <- rbind(deg_df_tmp, deg_sup_df)
}
unique(deg_sup_df$easyid_tumor)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Tumor_vs_PT_DEGs.", run_id, ".tsv")
write.table(x = deg_sup_df, file = file2write, quote = F, sep = "\t", row.names = F)
