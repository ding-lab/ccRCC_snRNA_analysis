# Yige Wu @WashU Sep 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(dplyr)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
gene2celltype_df <- fread(data.table = F, input = "./Resources/Knowledge/Kidney_Markers/Gene2CellType_Tab.20200911.v1.tsv")
emtgenes_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/write_emt_genes/20200911.v1/EMT_Genes.20200911.v1.tsv")

# combine -----------------------------------------------------------------

## prepare EMT genes
emtgenes_2combine_df <- emtgenes_df %>%
  dplyr::mutate(Gene_Group1 = "EMT_Genes") %>%
  dplyr::filter(Key_EMT_Genes) %>%
  # dplyr::filter(is.na(gene_group) | (!is.na(gene_group) & gene_group != "Collagen")) %>%
  dplyr::rename(Gene = hgnc_symbol) %>%
  dplyr::rename(Gene_Group2 = gene_function) %>%
  dplyr::select(Gene, Gene_Group1, Gene_Group2)
## prepare cell type markers
gene2celltype_2combine_df <- gene2celltype_df %>%
  dplyr::filter(!(Gene %in% emtgenes_2combine_df$Gene)) %>%
  dplyr::filter(Cell_Type1 %in% c("Tumor cells", "Proximal tubule", "Fibroblasts") ) %>%
  dplyr::mutate(Gene_Group1 = "Cell_Type_Markers") %>%
  dplyr::rename(Gene_Group2 = Cell_Type1) %>%
  dplyr::select(Gene, Gene_Group1, Gene_Group2)
## combine
combinedgenes_df <- rbind(emtgenes_2combine_df, gene2celltype_2combine_df)
  
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Kidney_Specific_EMT_Genes.", run_id, ".tsv")
write.table(x = combinedgenes_df, file = file2write, quote = F, sep = "\t", row.names = F)

