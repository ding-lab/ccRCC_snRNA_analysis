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
paths_monocle_objs <- fread(data.table = F, input = "./Resources/Analysis_Results/monocle/write_paths_to_monocle_objects/20200728.v1/Paths_to_Monocle_Objects.20200728.v1.tsv")

# process by case ---------------------------------------------------------
# for (id_run in paths_monocle_objs$Id_Run) {
for (id_run in "3Sample_PT_and_Tumorcells_ByCellTypeSampleID") {
  path_obj <- paths_monocle_objs$Path_Box[paths_monocle_objs$Id_Run == id_run]
  ## input monocle object
  obj_monocle <- readRDS(file = path_obj)
  ## plot
  p <- plot_cell_trajectory(obj_monocle, color_by = "Cell_type.detailed",cell_size=0.3)
  # p <- p + ggtitle(paste0("Pseudotime Trajectory of the Cells from\n", id_run))
  # p <- p + theme(aspect.ratio=1)
  p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
  ## write output
  file2write <- paste0(dir_out, id_run, ".png")
  png(file2write, width = 1000, height = 1200, res = 150)
  print(p)
  dev.off()
  file2write <- paste0(dir_out, id_run, ".pdf")
  pdf(file2write, width = 6, height = 6, useDingbats = F)
  print(p)
  dev.off()
}

