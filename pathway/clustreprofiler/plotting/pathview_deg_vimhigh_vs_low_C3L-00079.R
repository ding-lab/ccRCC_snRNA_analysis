# Yige Wu @WashU Sep 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(pathview)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
foldchanges_obj <- readRDS(file = "./Resources/Analysis_Results/pathway/clustreprofiler/clusterprofiler_deg_tumorlike_vs_tumorcells_C3L-00079/20200915.v1/Fold_Changes.RDS")

# plot --------------------------------------------------------------------
setwd(dir_out)
for (id_keggpathway in c("hsa04151")) {
  pathview(gene.data  = foldchanges_obj,
           pathway.id = id_keggpathway,
           species    = "hsa",
           limit      = list(gene=2, cpd=1))
}


