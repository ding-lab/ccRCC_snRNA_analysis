# Yige Wu @WashU Aug 2020

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
path_rds <- "./Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/subset_C3N-01200_tumorlikecells_and_recluster/20200910.v1/TumorLikeCells.Reclustered.20200910.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")

# specify genes to plot ---------------------------------------------------
# genes_plot <- c("VHL", "VIM", "COL5A1", "MT2A", "CDH1", "POSTN", "UMOD")
genes_plot <- c("VIM")

# plot by gene ------------------------------------------------------------
for (gene_plot in genes_plot) {
  p <- FeaturePlot(object = srat, features = gene_plot, order = T, min.cutoff = "q1", max.cutoff = "q99", cols = c("yellow", "red"))
  ## write output
  file2write <- paste0(dir_out, gene_plot, ".png")
  png(file2write, width = 1000, height = 1000, res = 150)
  print(p)
  dev.off()

}

