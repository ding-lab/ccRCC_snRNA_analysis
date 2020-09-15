# Yige Wu @WashU Sep 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(clusterProfiler)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
ora_obj <- readRDS(file = "./Resources/Analysis_Results/pathway/clustreprofiler/clusterprofiler_deg_tumorlike_vs_tumorcells_C3L-00079/20200910.v1/ORA_Results.RDS")

# plot --------------------------------------------------------------------
p <- dotplot(object = ora_obj)
p

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "dotplot.png")
png(file2write, width = 1200, height = 800, res = 150)
print(p)
dev.off()

