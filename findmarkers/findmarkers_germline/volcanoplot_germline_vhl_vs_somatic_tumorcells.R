# Yige Wu @WashU Jun 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(ggridges)
library(viridis)
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input denpendencies -----------------------------------------------------
## input marker table
markers_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_germline/findallmarker_wilcox_germline_vhl_vs_somatic_bycelltype_on_katmai/20200604.v1/VHL_Germline_vs_VHL_Somatic.FindMarkers.Wilcox.20200604.v1.tsv")
## input ppi table
ppi_pair_df <- fread(data.table = F, input = "./Resources/Databases/Protein_Protein_Interactions/protein_pair_table_v2.txt")

# get genes to highlight ------------------------------------------------------------------
## filter vhl interactome
genes_highlight_df <- ppi_pair_df %>%
  filter(GENE == "VHL")

# make plot data  ---------------------------------------------------------
celltype_plot <- "Tumor cells"
plot_data_df <- markers_df %>%
  filter(Cell_type.shorter == celltype_plot) %>%
  filter(p_val_adj < 0.05) %>%
  mutate(is.highlight = (deg_gene_symbol %in% genes_highlight_df$SUB_GENE))

# plot --------------------------------------------------------------------
p <- ggplot()
p <- p + geom_point(data = plot_data_df, mapping = aes(x = ))

# save output -------------------------------------------------------------

