# Yige Wu @WashU Mar 2022

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "plyr",
  "dplyr",
  "data.table"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    install.packages(pkg_name_tmp, dependencies = T)
  }
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg_name_tmp)
  }
  library(package = pkg_name_tmp, character.only = T)
}

# input dependencies ------------------------------------------------------
## input cell type per barcode table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups_35aliquots/20210802.v1/35Aliquot.Barcode2CellType.20210802.v1.tsv", data.table = F)
## input UMAP info per barcode
barcodes_keep_df <- fread(input = "./Resources/Analysis_Results/fetch_data/fetchdata_34_ccRCC_samples_merged_katmai/20211005.v1/ccRCC.34Sample.Merged.Metadata.20211005.v1.tsv", data.table = F)

# process -----------------------------------------------------------------
barcodes_info_df <- merge(barcodes_keep_df %>%
                        mutate(individual_barcode = str_split_fixed(string = barcode, pattern = "_", n = 2)[,1]),
                      barcode2celltype_df %>%
                        select(orig.ident, individual_barcode, Cell_type.shorter, Cell_type.detailed,
                               Cell_group_w_epithelialcelltypes),
                      by = c("orig.ident", "individual_barcode"), all.x = T)

# write output ------------------------------------------------------------
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
## make file name
file2write <- paste0(dir_out, "UMAPCoordinate.ClusterID.CellType.ByBarcode.tsv")
write.table(x = barcodes_info_df, file = file2write, quote = F, sep = "\t", row.names = F)
