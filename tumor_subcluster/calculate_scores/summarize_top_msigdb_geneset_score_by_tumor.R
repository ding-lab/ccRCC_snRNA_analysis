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
scores_wide_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/calculate_msigdb_top_geneset_scores/20210419.v1/MSigDB.Hallmark.tsv")

# summarize by sample -----------------------------------------------------
scores_long_df <- reshape2::melt(data = scores_wide_df, id.vars = c("cluster_name"))
scores_long_df <- scores_long_df %>%
  mutate(sample_id = str_split_fixed(string = cluster_name, pattern = "_", n = 2)[,1]) %>%
  group_by(sample_id, variable) %>%
  slice_max(order_by = value, n = 1)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Max.GenesetScorePerSample.tsv")
write.table(x = scores_long_df, file = file2write, quote = F, sep = "\t", row.names = F)
