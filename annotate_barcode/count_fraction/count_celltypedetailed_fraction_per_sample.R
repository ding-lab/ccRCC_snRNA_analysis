# Yige Wu @WashU Oct 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "data.table",
  "stringr"
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

# input data ------------------------------------------------
## input barcodes to cell type
barcode2celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups_35aliquots/20210802.v1/35Aliquot.Barcode2CellType.20210802.v1.tsv")
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210809.v1/meta_data.20210809.v1.tsv")

# specify Cell group to plot ----------------------------------------------
var_cellgroup <- "Cell_type.detailed"

# make data for plotting --------------------------------------------------
barcode2celltype_df$Aliquot_WU <- mapvalues(x = barcode2celltype_df$orig.ident, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
barcode2celltype_df[, "Cell_group"] <- barcode2celltype_df[, var_cellgroup]
## sum the barcode fraction by cell group
frac_barcodes_by_cellgroup <- barcode2celltype_df %>%
  group_by(Aliquot_WU, Cell_group) %>%
  summarise(Num_Barcodes_ByAliquotCellGroup = n())
num_barcodes <- barcode2celltype_df %>%
  group_by(Aliquot_WU) %>%
  summarise(Num_Barcodes_ByAliquot = n())
frac_barcodes_by_cellgroup$Num_Barcodes_ByAliquot <- mapvalues(x = frac_barcodes_by_cellgroup$Aliquot_WU, from = num_barcodes$Aliquot_WU, to = as.vector(num_barcodes$Num_Barcodes_ByAliquot))
frac_barcodes_by_cellgroup$Num_Barcodes_ByAliquot <- as.numeric(frac_barcodes_by_cellgroup$Num_Barcodes_ByAliquot)
frac_barcodes_by_cellgroup <- frac_barcodes_by_cellgroup %>%
  dplyr::mutate(Frac_CellGroupBarcodes_ByAliquot = Num_Barcodes_ByAliquotCellGroup/Num_Barcodes_ByAliquot)
## make data frame
plot_df <- frac_barcodes_by_cellgroup %>%
  mutate(Aliquot_Suffix = str_split_fixed(string = Aliquot_WU, pattern = "-", n = 3)[,3])
## label sample to tumor and normal
plot_df$Sample_Type <- mapvalues(x = plot_df$Aliquot_WU, from = idmetadata_df$Aliquot.snRNA.WU, to = idmetadata_df$Sample_Type)
plot_df$Case <- mapvalues(x = plot_df$Aliquot_WU, from = idmetadata_df$Aliquot.snRNA.WU, to = idmetadata_df$Case)

# write output ------------------------------------------------------------
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
file2write <- paste0(dir_out, "CellGroupBarcodes_Number_and_Fraction_per_Sample", run_id, '.tsv')
write.table(x = plot_df, file = file2write, quote = F, row.names = F, sep = "\t")
