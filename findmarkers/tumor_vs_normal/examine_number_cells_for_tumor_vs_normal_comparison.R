# Yige Wu @WashU Apr 2020
## summarize the number of cells will be used for tumor normal comparison

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
## input the barcode-cell-type table
path_barcode2celltype <- "./Resources/Analysis_Results/map_celltype_to_barcode/map_celltype_to_all_cells/20200410.v1/30_aliquot_integration.barcode2celltype.20200410.v1.tsv"
barcode2celltype_df <- fread(input = path_barcode2celltype, data.table = F)
## input id meta data
id_metadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20191105.v1/meta_data.20191105.v1.tsv", data.table = F)

# summarize the number of tumor cells and normal proximal tubule c --------
normalpt_barcode2celltype_df <- barcode2celltype_df %>%
  filter(Cell_type.detailed == "Proximal tubule")
normalpt_barcode2celltype_df$id_case <- mapvalues(x = normalpt_barcode2celltype_df$orig.ident, from = id_metadata_df$Aliquot.snRNA, to = as.vector(id_metadata_df$Case))
count_normalpt_byaliquot_df <- normalpt_barcode2celltype_df %>%
  select(orig.ident, id_case) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq > 0) %>%
  rename(id_aliquot = orig.ident)
count_normalpt_bycase_df <- normalpt_barcode2celltype_df %>%
  select(id_case) %>%
  table() %>%
  as.data.frame() %>%
  rename(id_case = '.')
count_tumorcells_df <- barcode2celltype_df %>%
  filter(Cell_type.detailed == "Tumor cells") %>%
  select(orig.ident) %>%
  table() %>%
  as.data.frame() %>%
  rename(id_aliquot = '.')
count_tumorcells_df$id_case <- mapvalues(x = count_tumorcells_df$id_aliquot, from = id_metadata_df$Aliquot.snRNA, to = as.vector(id_metadata_df$Case))

# merge tumor cell and normal epithelial cell number ----------------------
count_tumorcells_normalpt_bycase_df <- merge(count_tumorcells_df, count_normalpt_bycase_df, by = c("id_case"), all.y = T, suffixes = c(".tumorcells", ".normalpt"))
count_tumorcells_normalpt_bycase_df <- count_tumorcells_normalpt_bycase_df %>%
  arrange(-Freq.normalpt)

count_tumorcells_normalpt_byaliquot_df <- merge(count_tumorcells_df, count_normalpt_byaliquot_df, by = c("id_aliquot", "id_case"), all.y = T, suffixes = c(".tumorcells", ".normalpt"))
count_tumorcells_normalpt_byaliquot_df <- count_tumorcells_normalpt_byaliquot_df %>%
  filter(!is.na(Freq.tumorcells)) %>%
  arrange(-Freq.normalpt)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "number_of_malignant_paired_with_normal_proximal_tubule_cells_by_case.", run_id, ".tsv")
write.table(x = count_tumorcells_normalpt_bycase_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "number_of_malignant_paired_with_normal_proximal_tubule_cells_by_aliquot.", run_id, ".tsv")
write.table(x = count_tumorcells_normalpt_byaliquot_df, file = file2write, quote = F, sep = "\t", row.names = F)

