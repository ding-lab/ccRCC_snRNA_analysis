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
enrich_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/pathway/gsea_tumor_manualsubcluster_degs/20210413.v1/GSEA_Results.tsv")

# count -------------------------------------------------------------------
count_geneset_df <- enrich_df %>%
  filter(pvalue < 0.05) %>%
  dplyr::select(Description, easy_id) %>%
  unique() %>%
  dplyr::select(Description) %>%
  table() %>%
  as.data.frame() %>%
  dplyr::rename(GeneSet_Name = ".") %>%
  arrange(desc(Freq))
file2write <- paste0(dir_out, "Count.GSEA.Pvalue.0.05.tsv")
write.table(x = count_geneset_df, file = file2write, quote = F, sep = "\t", row.names = F)

count_geneset_df <- enrich_df %>%
  filter(p.adjust < 0.05) %>%
  dplyr::select(Description, easy_id) %>%
  unique() %>%
  dplyr::select(Description) %>%
  table() %>%
  as.data.frame() %>%
  dplyr::rename(GeneSet_Name = ".") %>%
  arrange(desc(Freq))
file2write <- paste0(dir_out, "Count.GSEA.P.adjust.0.05.tsv")
write.table(x = count_geneset_df, file = file2write, quote = F, sep = "\t", row.names = F)
