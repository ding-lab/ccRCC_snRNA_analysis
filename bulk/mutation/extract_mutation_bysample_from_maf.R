# Yige Wu @WashU Dec 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/load_data.R")

## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the bulk mutation data
maf_df <- loadMaf()

# preprocess --------------------------------------------------------------
maf_df <- maf_df %>%
  mutate(Case = str_split_fixed(string = Tumor_Sample_Barcode, pattern = "_", n = 2)[,1])
table(maf_df$Case)

# write output ------------------------------------------------------------
case_tmp <- "C3L-00096"
maf_case_df <- maf_df %>%
  filter(Case == case_tmp)

file2write <- paste0(dir_out, case_tmp, ".ccrcc.somatic.consensus.gdc.umichigan.wu.112918", ".tsv")
write.table(x = maf_case_df, file = file2write, quote = F, sep = "\t", row.names = F)
