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

# set parameters for findmarkers ------------------------------------------
## set seurat ids for group 1 and 2
clusterids_tumorlike <- 7
clusterids_tumorcells <- c(10, 4, 1, 0, 6)
## parameters for the findmarkers
# logfc.threshold.run <- log(2)
logfc.threshold.run <- 0
min.pct.run <- 0.1
min.diff.pct.run <- 0.1
## spcify assay
assay_process <- "RNA"
cat(paste0("Assay: ", assay_process, "\n"))
cat("###########################################\n")

# run findallmarkers ------------------------------------------------------
degs_df <- FindMarkers(object = srat, test.use = "wilcox", ident.1 = clusterids_tumorlike, ident.2 = clusterids_tumorcells, only.pos = F,
                             min.pct = min.pct.run, logfc.threshold = logfc.threshold.run, min.diff.pct = min.diff.pct.run, verbose = T)
degs_df$row_name <- rownames(degs_df)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "findmarkers_wilcox_tumorlikecells_vs_tumorcells", ".logfcthreshold", logfc.threshold.run, ".minpct", min.pct.run, ".mindiffpct", min.diff.pct.run, ".tsv")
write.table(x = degs_df, file = file2write, sep = "\t", quote = F, row.names = F)

