# Yige Wu @WashU Jul 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(monocle)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# input dependencies ------------------------------------------------------
## input monocle object
obj_monocle <- readRDS(file = "./Resources/snRNA_Processed_Data/Monocle/outputs/CellTypeVer.20200828.v1/C3L-00079_PooledPT_SelfFib_ByCellType/combined_subset_pseudotime_qval_1e-10.rds")

# specify genes to plot ---------------------------------------------------
genes_plot <- c("VIM", "CDH1")

# process by case ---------------------------------------------------------
gene_plot <- "VIM"
## 
p <- plot_genes_in_pseudotime(obj_monocle[gene_plot,], color_by = "Cell_type.detailed")
file2write <- paste0(dir_out, "plot_genes_in_pseudotime", ".png")
png(file2write, width = 1000, height = 1000, res = 150)
print(p)
dev.off()

plot_genes_branched_pseudotime(cds = obj_monocle[genes_plot,], branch_point = 1, color_by = "Cell_type.detailed", )
file2write <- paste0(dir_out, "plot_genes_branched_pseudotime", ".png")
png(file2write, width = 1000, height = 2000, res = 150)
print(p)
dev.off()


