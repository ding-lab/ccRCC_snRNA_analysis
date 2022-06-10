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

# input ------------------------------------------------------
## input the manual intrapatient cluster
barcode2tumorsubcluster_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_tumorsubclusterid/map_barcode_with_manual_tumorsubcluster_id/20210805.v1/Barcode2TumorSubclusterId.20210805.v1.tsv", data.table = F)
## input the tumor cell intergrated barcode info
barcode_metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/integration/seuratintegrate_34_ccRCC_samples/FindClusters_30_ccRCC_tumorcells_changeresolutions/20220405.v1/ccRCC.34Sample.Tumorcells.Integrated.ReciprocalPCA.Metadata.ByResolution.20220405.v1.tsv")

# pre-process -----------------------------------------------------------------
barcode2tumorsubcluster_df <- barcode2tumorsubcluster_df %>%
  mutate(barcode_id_uniq = paste0(orig.ident, "_", barcode))
barcode_metadata_df <- barcode_metadata_df %>%
  mutate(barcode_individual = str_split_fixed(string = barcode, pattern = "_", n = 2)[,1]) %>%
  mutate(barcode_id_uniq = paste0(orig.ident, "_", barcode_individual))
barcode_metadata_df <- merge(x = barcode_metadata_df, 
                             y = barcode2tumorsubcluster_df %>%
                               select(barcode_id_uniq, Cluster_Name) %>%
                               rename(intrapatient_cluster_name = Cluster_Name), 
                             by = c("barcode_id_uniq"), all.x = T)
barcode_bymc_byic_df <- barcode_metadata_df %>%
  mutate(meta_cluster_name = paste0("MC", integrated_snn_res.1)) %>%
  group_by(intrapatient_cluster_name, meta_cluster_name) %>%
  summarise(cells_bymc_byic = n())
barcode_byic_df <- barcode_metadata_df %>%
  group_by(intrapatient_cluster_name) %>%
  summarise(cells_byic = n())
barcode_bymc_byic_df <- merge(x = barcode_bymc_byic_df, y = barcode_byic_df, by = c("intrapatient_cluster_name"), all.x = T)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "median_vision_score.sig_genesets.byintrapatient_tumor_cluster.", run_id,".tsv")
write.table(x = barcode_bymc_byic_df, file = file2write, quote = F, sep = "\t", row.names = F)


