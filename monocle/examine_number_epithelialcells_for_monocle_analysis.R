# Yige Wu @WashU Apr 2020
## summarize the number of cells will be used for tumor normal comparison

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

# input dependencies ------------------------------------------------------
## input the barcode-cell-type table
path_barcode2celltype <- "./Resources/Analysis_Results/annotate_barcode/map_celltype_corrected_by_individual_sample_inspection/20200811.v1/31Aliquot.Barcode2CellType.20200811.v1.tsv"
barcode2celltype_df <- fread(input = path_barcode2celltype, data.table = F)
## input id meta data
id_metadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv", data.table = F)

# map case id -------------------------------------------------------------
barcode2celltype_df$id_case <- mapvalues(x = barcode2celltype_df$orig.ident, from = id_metadata_df$Aliquot.snRNA, to = as.vector(id_metadata_df$Case))

# count by cell type shorter ---------------------------------------------
count_epicelltypeshorter_bycase_long_df <- barcode2celltype_df %>%
  dplyr::filter(Cell_group.shorter == "Nephron_Epithelium") %>%
  dplyr::select(id_case, Cell_type.shorter) %>%
  table() %>%
  data.frame()
count_epicelltypeshorter_bycase_wide_df <- dcast(data = count_epicelltypeshorter_bycase_long_df, formula = id_case~Cell_type.shorter, value.var = "Freq")
# count by cell type detailed ---------------------------------------------
count_epicelltypedetailed_bycase_long_df <- barcode2celltype_df %>%
  dplyr::filter(Cell_group.shorter == "Nephron_Epithelium") %>%
  dplyr::select(id_case, Cell_type.detailed) %>%
  table() %>%
  data.frame()
count_epicelltypedetailed_bycase_wide_df <- dcast(data = count_epicelltypedetailed_bycase_long_df, formula = id_case~Cell_type.detailed, value.var = "Freq")

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "number_of_epithelialcells_bycelltypeshorter_by_case.", run_id, ".tsv")
write.table(x = count_epicelltypeshorter_bycase_wide_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "number_of_epithelialcells_bycelltypedetailed_by_case.", run_id, ".tsv")
write.table(x = count_epicelltypedetailed_bycase_wide_df, file = file2write, quote = F, sep = "\t", row.names = F)

