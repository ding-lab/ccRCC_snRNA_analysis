# Yige Wu @WashU Sep 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(clusterProfiler)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
ora_obj <- readRDS(file = "./Resources/Analysis_Results/pathway/clustreprofiler/clusterprofiler_deg_tumorlike_vs_tumorcells_C3L-00079/20200915.v1/ORA_Results.RDS")
foldchanges_obj <- readRDS(file = "./Resources/Analysis_Results/pathway/clustreprofiler/clusterprofiler_deg_tumorlike_vs_tumorcells_C3L-00079/20200915.v1/Fold_Changes.RDS")

# plot --------------------------------------------------------------------
p <- cnetplot(ora_obj, foldChange=foldchanges_obj, showCategory = 3)
p <- p + scale_color_gradient2(low = "blue", high = "red")
p <- p + labs(colour = "Fold change\n(transitional vs\ntumor cells)")
p

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "gene_concept_netowrk.png")
png(file2write, width = 1200, height = 1000, res = 150)
print(p)
dev.off()

file2write <- paste0(dir_out, "gene_concept_netowrk.pdf")
pdf(file2write, width = 8, height = 6, useDingbats = F)
print(p)
dev.off()

