# Yige Wu @WashU Apr 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
enrich_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/pathway/ora_msigdb_tumor_manualsubcluster_up_degs/20210413.v1/ORA_Results.tsv")
enrich_df <- fread(data.table = F, input = "~/Downloads/ORA_Results.tsv")

# make plot data ----------------------------------------------------------
count_geneset_df <- enrich_df %>%
  filter(p.adjust < 0.05) %>%
  dplyr::select(Description, easy_id) %>%
  unique() %>%
  dplyr::select(Description) %>%
  table() %>%
  as.data.frame() %>%
  dplyr::rename(Description = ".") %>%
  arrange(desc(Freq))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Count_gene_set_in_up_tumorcluster_degs.", run_id, ".tsv")
write.table(x = count_geneset_df, file = file2write, quote = F, sep = "\t", row.names = F)
