# Yige Wu @WashU Mar 2022

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "Seurat",
  "ggrastr",
  "ggplot2",
  "readxl"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
source("./ccRCC_snRNA_analysis/variables.R")

# input dependencies ------------------------------------------------------
## input 10xmapping result
snRNA_mutation_df <- fread("./Resources/Analysis_Results/mutation/unite_10xmapping/20220126.v1/10XMapping.20220126.v1.tsv", data.table = F)
barcode2tumorsubcluster_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_tumorsubclusterid/map_barcode_with_manual_tumorsubcluster_id/20210805.v1/Barcode2TumorSubclusterId.20210805.v1.tsv", data.table = F)
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20210729.v1/scrublet.united_outputs.20210729.v1.tsv", data.table = F)

# pre-process -------------------------------------------------------------
barcode2scrublet_df <- barcode2scrublet_df %>%
  mutate(barcode_uniq = paste0(Aliquot, "_", Barcode))
barcodes_doublet <- barcode2scrublet_df$barcode_uniq[barcode2scrublet_df$predicted_doublet]
driver_mutation_df <- snRNA_mutation_df %>%
  filter(gene_symbol %in% ccRCC_drivers) %>%
  filter(allele_type == "Var") %>%
  mutate(barcode_uniq = paste0(aliquot, "_", barcode))
barcodes_mapped <- driver_mutation_df$barcode_uniq
barcode2tumorsubcluster_df <- barcode2tumorsubcluster_df%>%
  mutate(barcode_uniq = paste0(orig.ident, "_", barcode)) %>%
  filter(!(barcode_uniq %in% barcodes_doublet)) %>%
  mutate(driver_mutation_mapped = (barcode_uniq %in% barcodes_mapped))
nrow(barcode2tumorsubcluster_df)
driver_mutation_bytumorcluster_df <- barcode2tumorsubcluster_df %>%
  group_by(Cluster_Name) %>%
  summarise(number_cells_w_driver_mutation = length(which(driver_mutation_mapped))) %>%
  arrange(desc(number_cells_w_driver_mutation))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Driver_mutation_mapped_per_intrapatienttumorcluster.", run_id, ".tsv")
write.table(x = driver_mutation_bytumorcluster_df, file = file2write, quote = F, sep = "\t", row.names = F)
