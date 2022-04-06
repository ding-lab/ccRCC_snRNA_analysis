# Yige Wu @WashU Apr 2022
## plot cell type on integration UMAP

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input UMAP info per barcode
barcode2umap_df <- fread(input = "./Resources/Analysis_Results/integration/seuratintegrate_34_ccRCC_samples/reciprocalPCA_integrate_30_ccRCC_tumorcells/20220404.v1/ccRCC.34Sample.Tumorcells.Integrated.ReciprocalPCA.Metadata.20220404.v1.tsv", data.table = F)
## input barcode to cluster id
barcode2cluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/integration/seuratintegrate_34_ccRCC_samples/FindClusters_30_ccRCC_tumorcells_changeresolutions/20220405.v1/ccRCC.34Sample.Tumorcells.Integrated.ReciprocalPCA.Metadata.ByResolution.20220405.v1.tsv")

# make plot data----------------------------------------------------------
barcode_info_df <- merge(x = barcode2umap_df %>%
                        select(barcode, UMAP_1, UMAP_2),
                      y = barcode2cluster_df %>%
                        select(orig.ident, barcode, integrated_snn_res.0.1, integrated_snn_res.0.2, integrated_snn_res.0.3, integrated_snn_res.0.4, integrated_snn_res.0.5,
                               integrated_snn_res.1, integrated_snn_res.2),
                      by = c("barcode"), all.x = T)
table(barcode_info_df$integrated_snn_res.0.5)
## Merge cluster 0, 1, 2, 4, 6 together, while keeping the rest separate
barcode_info_df$clusterid_new <- mapvalues(x = barcode_info_df$integrated_snn_res.0.5, 
                                           from = c(0, 1, 2, 4, 6,
                                                    3, 5, 7, 8, 9, 10, 11), 
                                           to = c("C1", "C1", "C1", "C1", "C1", 
                                                  "C2", "C3", "C4", "C5", "C6", "C7", "C8"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "ccRCC.34Sample.Tumorcells.Integrated.ReciprocalPCA.Metadata.", run_id, ".tsv")
write.table(x = barcode_info_df, file = file2write, sep = "\t", row.names = F, quote = F)


