# Yige Wu @WashU Feb 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
srat <- readRDS(file = "./Resources/Analysis_Results/integration/merge_different_cases/merge_tumorcells_C3L-00010_C3L-00416_C3L-00583_C3N-00242_on_katmai/20200817.v1/C3L-00416_C3L-00583_C3L-00010_C3N-00242.Tumor_Segments.Merged.20200817.v1.RDS")
print("Finish reading the RDS file!")

# set parameters for findmarkers ------------------------------------------
## set min.pct
logfc.threshold.wilcox <- 0.1
min.pct.wilcox <- 0.1

# find all markers --------------------------------------------------------
markers_wilcox_df <- FindAllMarkers(object = srat, test.use = "wilcox", only.pos = F, min.pct = min.pct.wilcox, logfc.threshold = logfc.threshold.wilcox, verbose = T)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "C3N-00010.correlatedtumorcells.FindAllMarkers.bynewcluster.Wilcox.Minpct", min.pct.wilcox, ".Logfc", logfc.threshold.wilcox,".tsv")
write.table(markers_wilcox_df, file = file2write, 
            quote = F, sep = "\t", row.names = F)