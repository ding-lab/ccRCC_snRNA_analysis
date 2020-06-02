# Yige Wu @WashU May 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(ggpubr)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies -------------------------------------------------------
## input cnv plot
cnv_plot <- readRDS(file = "./Resources/Analysis_Results/copy_number/summarize_cnv_fraction/bubble_cnv_fraction_per_tumormanualsubcluster/20200602.v2/bubbleplot_cnv_fraction.20200602.v2.RDS")
## input deg plot
deg_plot <- readRDS(file = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/bubbleplot_deg_by_tumormanualsubcluster/20200602.v1/bubbleplot_deg_by_tumorsubcluster.20200602.v1.RDS")

# combine plots together --------------------------------------------------
# p <- ggarrange(cnv_plot, deg_plot, align = "v", 
#                labels = c("A", "B"),
#                ncol = 1, nrow = 2)
p <- cowplot::plot_grid(cnv_plot, deg_plot, ncol = 1, align = "v", axis = "l", rel_heights = c(3,2))


# save output -------------------------------------------------------------
## save plot
file2write <- paste0(dir_out, "bubbleplot_deg_by_tumorsubcluster.", run_id, ".png")
png(file = file2write, 
    width = 2500, height = 2000, res = 150)
print(p)
dev.off()


