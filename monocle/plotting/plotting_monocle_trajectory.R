# Yige Wu @WashU Jul 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input monocle object
obj_monocle <- readRDS(file = "./Resources/snRNA_Processed_Data/Monocle/outputs/C3N-01200/combined_subset_pseudotime_qval_1e-10.rds")
## input cell type and tumor subcluster info
barcode2celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_to_all_cells_with_patch/20200720.v1/31AliquotIntegration.Barcode2CellType.TumorManualCluster.20200720.v1.tsv")


# plot --------------------------------------------------------------------

plot_cell_trajectory(cds = obj_monocle, color_by = 'Cell_type.shorter')


plot_cell_trajectory(cds = obj_monocle, color_by = 'Naming')
