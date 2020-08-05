# Yige Wu @WashU Jul 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(monocle)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input cell type
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_to_all_cells/20200720.v1/30AliquotIntegration.Barcode2CellType.TumorManualCluster.20200720.v1.tsv", data.table = F)
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200505.v1/meta_data.20200505.v1.tsv")

# add case id -------------------------------------------------------------
barcode2celltype_df$Id_Case <- mapvalues(x = barcode2celltype_df$orig.ident, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Case))

# check C3N-01200 ---------------------------------------------------------
barcode2celltype_df %>%
  filter(Id_Case == "C3N-01200") %>%
  filter(Most_Enriched_Cell_Group == "Nephron_Epithelium") %>%
  select(Cell_type.detailed) %>%
  table()
barcode2celltype_df %>%
  filter(Id_Case == "C3N-01200") %>%
  filter(Most_Enriched_Cell_Group == "Nephron_Epithelium") %>%
  select(Cell_group) %>%
  table()
## I'll do one run with all epithelial cells and on with 1000 tumor cells



