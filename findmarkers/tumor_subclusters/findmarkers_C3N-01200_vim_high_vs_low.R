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
path_rds <- "./Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/subset_C3N-01200_tumorlikecells_and_recluster/20200910.v1/TumorLikeCells.Reclustered.20200910.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")

# set parameters for findmarkers ------------------------------------------
## set seurat ids for group 1 and 2
clusterids_group1 <- c(1, 6)
clusterids_group2 <- c(0, 2, 3, 4, 5)
## parameters for the findmarkers
logfc.threshold.run <- log(2)
min.pct.run <- 0.1
min.diff.pct.run <- 0.1
## spcify assay
assay_process <- "RNA"
cat(paste0("Assay: ", assay_process, "\n"))
cat("###########################################\n")

# run findallmarkers ------------------------------------------------------
degs_df <- NULL
for (clusterid_group1 in clusterids_group1) {
  degs_tmp_df <- FindMarkers(object = srat, test.use = "wilcox", ident.1 = clusterid_group1, ident.2 = clusterids_group2, only.pos = F,
                         min.pct = min.pct.run, logfc.threshold = logfc.threshold.run, min.diff.pct = min.diff.pct.run, verbose = T)
  degs_tmp_df$row_name <- rownames(degs_tmp_df)
  degs_tmp_df$clusterid_group1 <- clusterid_group1
  degs_df <- rbind(degs_df, degs_tmp_df)
}

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "findmarkers_wilcox_vimhigh_vs_low", ".logfcthreshold", logfc.threshold.run, ".minpct", min.pct.run, ".mindiffpct", min.diff.pct.run, ".tsv")
write.table(x = degs_df, file = file2write, sep = "\t", quote = F, row.names = F)

