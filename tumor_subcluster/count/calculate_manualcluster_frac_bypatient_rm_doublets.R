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
  "data.table"
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

# input dependencies ------------------------------------------------------
barcode2tumorsubcluster_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_tumorsubclusterid/map_barcode_with_manual_tumorsubcluster_id/20210805.v1/Barcode2TumorSubclusterId.20210805.v1.tsv", data.table = F)
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20210729.v1/scrublet.united_outputs.20210729.v1.tsv", data.table = F)

# pre-process --------------------------------------------------------------
barcode2tumorsubcluster_df <- barcode2tumorsubcluster_df %>%
  mutate(case = str_split_fixed(string = easy_id, pattern = "\\-T", n = 2)[,1]) %>%
  mutate(barcode_uniq = paste0(orig.ident, "_", barcode))
barcode2scrublet_df <- barcode2scrublet_df %>%
  mutate(barcode_uniq = paste0(Aliquot, "_", Barcode))
barcodes_doublet <- barcode2scrublet_df$barcode_uniq[barcode2scrublet_df$predicted_doublet]
barcode2tumorsubcluster_filtered_df <- barcode2tumorsubcluster_df %>%
  filter(!(barcode_uniq %in% barcodes_doublet))

# count ----------------------------------------------------
count_bycluster_df <- barcode2tumorsubcluster_filtered_df %>%
  group_by(Cluster_Name, easy_id, case) %>%
  summarise(count_bycluster = n())
count_bycase_df  <- barcode2tumorsubcluster_filtered_df %>%
  group_by(case) %>%
  summarize(count_bycase = n())
count_bycluster_df <- merge(x = count_bycluster_df, y = count_bycase_df, by = c("case"), all.x = T)
count_bycluster_df <- count_bycluster_df %>%
  mutate(frac_cluster_count_bycase = count_bycluster/count_bycase)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Fraction_cluster_by_case.", run_id, ".tsv")
write.table(x = count_bycluster_df, file = file2write, quote = F, sep = "\t", row.names = F)
