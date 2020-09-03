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
obj_monocle <- readRDS(file = "./Resources/snRNA_Processed_Data/Monocle/outputs/CellTypeVer.20200825.v1/C3L-00088_SelfEpiCells_ByCellType/combined_subset_pseudotime_qval_1e-10.rds")

# process by case ---------------------------------------------------------
## make colors
pdata_df <- data.frame(obj_monocle@phenoData@data)
aliquot_ids <- unique(pdata_df$Aliquot.snRNA.WU)
colors_aliquot <- c("orange", "purple", rep("grey50", (length(aliquot_ids) - 2)))
names(colors_aliquot) <- c("C3L-00088-T1", "C3L-00088-T2", as.vector(aliquot_ids[!(aliquot_ids %in% c("C3L-00088-T1", "C3L-00088-T2"))]))
## plot
p <- plot_cell_trajectory(obj_monocle, color_by = "Aliquot.snRNA.WU",cell_size=0.3)
p <- p + scale_color_manual(values = colors_aliquot)
# p <- p + ggtitle(paste0("Pseudotime Trajectory of the Cells from\n", id_run))
# p <- p + theme(aspect.ratio=1)
p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
## write output
file2write <- paste0(dir_out, "celltypedetailed", ".png")
png(file2write, width = 1000, height = 1000, res = 150)
print(p)
dev.off()
file2write <- paste0(dir_out, "Aliquot.snRNA.WU", ".pdf")
pdf(file2write, width = 6, height = 6, useDingbats = F)
print(p)
dev.off()

