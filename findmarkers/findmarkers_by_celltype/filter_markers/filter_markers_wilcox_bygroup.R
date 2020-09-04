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
## input the marker genes
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_by_celltype/findallmarker_wilcox_cellgroup_on_katmai/20200714.v2/findallmarkers_wilcox_bycellgroup.pos..logfcthreshold0.1.minpct0.1.mindiffpct0.1.tsv")

# get top 50 genes by avg_logFC ------------------------------------------------
unique(deg_df$cluster)
deg_filtered_df <- deg_df %>%
  group_by(cluster) %>%
  top_n(wt = avg_logFC, n = 50)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "findallmarkers_wilcox_bycellgroup.pos.logfcthreshold0.1.minpct0.1.mindiffpct0.1.top50avg_logFC",".tsv")
write.table(x = deg_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)

