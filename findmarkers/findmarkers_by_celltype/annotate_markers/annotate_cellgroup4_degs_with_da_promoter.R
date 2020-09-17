# Yige Wu @WashU Aug 2020

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

# input dependencies -----------------------------------------------------
genes_df <- fread(input = "./Resources/Analysis_Results/findmarkers/findmarkers_by_celltype/filter_markers/filter_markers_wilcox_bygroup/20200916.v1/findallmarkers_wilcox_bycellgroup.pos.logfcthreshold0.1.minpct0.1.mindiffpct0.1.Top50avg_logFC.tsv", data.table = F)
## input tumor vs normal DA peaks
peaks2celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/filter_peaks/filter_celltype4_da_peaks_for_cellgoup4/20200916.v1/Cell_Type_Specific_DEGs_in_DA_peaks.20200916.v1.tsv")
## specify top n for the deg table
n_top <- 50

# annotate gene to different levels of evidence ---------------------------
genes_df$CellTypeSpecific_DA_Promoter <- (genes_df$gene %in% peaks2celltype_df$SYMBOL)
genes_df <- genes_df %>%
  group_by(cluster) %>%
  mutate(rank_deg_by_avglogFC = order(order(avg_logFC, decreasing = T)))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "CellTypeDEGTop", n_top,  ".Annotation.", run_id, ".tsv")
write.table(x = genes_df, file = file2write, quote = F, sep = "\t", row.names = F)

