# Yige Wu @WashU Sep 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the integrated data
path_rds <- "./Resources/snRNA_Processed_Data/scRNA_auto/outputs/CPT0001260013/pf1000_fmin200_fmax10000_cmin1000_cmax10000_mito_max0.1/CPT0001260013_processed.rds"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")

# Dimplot -----------------------------------------------------------------
## make distinguishable colors
colors_cluster <- Polychrome::dark.colors(n = 15)
names(colors_cluster) <- 0:14
p <- DimPlot(object = srat)
p <- p + scale_color_manual(values = colors_cluster)
p
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Seurat_Cluster_Id.", "C3L-00079", ".png")
png(file2write, width = 900, height = 800, res = 150)
print(p)
dev.off()
