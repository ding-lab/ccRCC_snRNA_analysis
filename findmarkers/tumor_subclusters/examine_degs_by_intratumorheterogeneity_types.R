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

# input denpendencies -----------------------------------------------------
## input log fold change
min.pct.wilcox <- 0.1
logfc.threshold.wilcox <- 0.25
logfc_by_manualsubcluster_df <- fread(input = paste0("./Resources/Analysis_Results/findmarkers/tumor_subclusters/findallmarkers_wilcox_tumorcells_by_manualsubcluster/20200427.v1/Tumormanualsubcluster.FindAllMarkers.Wilcox.Minpct", min.pct.wilcox, ".Logfc", logfc.threshold.wilcox, ".tsv"), data.table = F)
## input phenotype markers
markers_ith_df <- fread("./Resources/Analysis_Results/dependencies/write_markers_by_intratumorheterogeneity_types/20200427.v1/markergenes_by_intratumorheterogeneity_types.20200427.v1.tsv", data.table = F)
## input correlation with CNV
spearman_result_df <- fread(input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/calc_cor_btw_deg_avgexp_and_frac_cnv/20200428.v1/spearman_cor_btw_deg_avgexp_and_frac_cnv.20200428.v1.tsv", data.table = F)

# intersect deg with phenotype markers---------------------------------------------------------------
logfc_by_manualsubcluster_df <- merge(logfc_by_manualsubcluster_df, markers_ith_df, by.x = c("gene"), by.y = c("gene_symbol"), all.x = T)
logfc_by_manualsubcluster_df <- logfc_by_manualsubcluster_df %>%
  filter(p_val_adj < 0.05) %>%
  filter(avg_logFC > 0)

# interesect correlation result gene with phenotype markers --------------------------
spearman_result_sig_df <- spearman_result_df %>%
  filter(spearman_fdr < 0.1)
## merge
spearman_result_sig_df <- merge(spearman_result_sig_df, markers_ith_df, by.x = c("gene_symbol"), by.y = c("gene_symbol"), all.x = T)

