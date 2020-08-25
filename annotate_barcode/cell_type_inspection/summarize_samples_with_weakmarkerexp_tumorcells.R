# Yige Wu @WashU Aug 2020

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
## input cell type per barcode table
barcode2celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_corrected_by_individual_sample_inspection/20200824.v1/31Aliquot.Barcode2CellType.20200824.v1.tsv")
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")

# count the interested tumor cells by sample ------------------------------
## map readable ids
barcode2celltype_df$id_aliquot_wu <- mapvalues(x = barcode2celltype_df$orig.ident, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
## count by case
barcode2celltype_filtered_df <- barcode2celltype_df %>%
  filter(Cell_type.detailed == "Tumor cells (weak marker gene expression)") %>%
  select(id_aliquot_wu) %>%
  table() %>%
  as.data.frame() %>%
  rename(Aliquot = '.') %>%
  arrange(desc(Freq))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Count_Tumorcells_w.weak_marker_gene_exp.", run_id, ".tsv")
write.table(x = barcode2celltype_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)

