# Yige Wu @WashU June 2020

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
## input degs
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/findallmarkers_wilcox_tumorcells_by_manualsubcluster_with_cnv_diff_on_katmai/20200623.v2/Tumormanualsubcluster.withCNVDiff.FindMarkers.Wilcox.Minpct0.Logfc0.1.tsv")
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200505.v1/meta_data.20200505.v1.tsv")
## input ppi table
ppi_pair_df <- fread(data.table = F, input = "./Resources/Knowledge/Databases/Protein_Protein_Interactions/protein_pair_table_v2.txt")
## set sample id
id_aliquot_wu_filter <- "C3L-00004-T1"
# filter by sample --------------------------------------------------------
deg_sample_df <- deg_df %>%
  filter(id_aliquot_wu == id_aliquot_wu_filter) %>%
  unique()
# filter ppit -------------------------------------------------------------
ppi_gene_filtered_df <- ppi_pair_df %>%
  filter(GENE %in% c("HIF1A", "ARID1A"))

# annotate degs -----------------------------------------------------------
ppi_pair_filtered_df <- ppi_gene_filtered_df %>%
  filter(SUB_GENE %in% degs_sig)
ppi_pair_filtered_df$SUB_GENE
deg_sample_filtered_df <- merge(deg_sample_df %>%
                                  filter(p_val < 0.05), ppi_pair_filtered_df, by.x = c("deg_gene_symbol"), by.y = c("SUB_GENE"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, id_aliquot_wu_filter, ".DEG_selected.tsv")
write.table(x = deg_sample_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)


