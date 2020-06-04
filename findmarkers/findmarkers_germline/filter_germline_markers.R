# Yige Wu @WashU Jun 2020

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

# input denpendencies -----------------------------------------------------
## input ppi table
ppi_pair_df <- fread(data.table = F, input = "./Resources/Databases/Protein_Protein_Interactions/protein_pair_table_v2.txt")
## input marker table
markers_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_germline/findallmarker_wilcox_germline_vhl_vs_somatic_tumorcells_on_katmai/20200604.v1/Germline_vs_Somatic_VHL.TumorCells.FindMarkers.Wilcox.20200604.v1.tsv")
markers_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_germline/findallmarker_wilcox_germline_vhl_vs_somatic_bycelltype_on_katmai/20200604.v1/VHL_Germline_vs_group2_findmarkers.FindMarkers.Wilcox.20200604.v1.tsv")

# filter markers ----------------------------------------------------------
markers_filtered_df <- markers_df %>%
  filter(p_val_adj < 0.05)

# filter vhl interactome ------------------------------------------------------------------
ppi_pair_filtered_df <- ppi_pair_df %>%
  filter(GENE == "VHL")

# merge -------------------------------------------------------------------
markers_filtered_anno_df <- merge(markers_filtered_df, ppi_pair_filtered_df, by.x = c("deg_gene_symbol"), by.y = c("SUB_GENE"), all.x = T)
markers_filtered_anno_df <- markers_filtered_anno_df %>%
  arrange(GENE, SUB_GENE.is_complex_partner, avg_logFC)

