# Yige Wu @WashU June 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input degs
# deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/findallmarkers_wilcox_tumorcells_by_manualsubcluster_with_cnv_diff_on_katmai/20200623.v2/Tumormanualsubcluster.withCNVDiff.FindMarkers.Wilcox.Minpct0.Logfc0.1.tsv")
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/findallmarkers_wilcox_tumorcells_by_manualsubcluster_with_cnv_maxdiff/20200603.v1/Tumormanualsubcluster.FindMarkers.Cluster_CNVLow_vs_CNVHigh.Wilcox.Minpct0.1.Logfc0.1.tsv")
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200505.v1/meta_data.20200505.v1.tsv")
## input ppi table
ppi_pair_df <- fread(data.table = F, input = "./Resources/Knowledge/Databases/Protein_Protein_Interactions/protein_pair_table_v2.txt")
## set sample id
id_aliquot_wu_filter <- "C3N-00242-T1"
# filter by sample --------------------------------------------------------
deg_sample_df <- deg_df %>%
  filter(id_aliquot_wu == id_aliquot_wu_filter) %>%
  select(-cna_gene_symbol) %>%
  unique()
# "CECR2" is the most significant DEG,a member of CERF chromatin-remodelling complex, which is important for CNS development: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3060774/
# filter ppit -------------------------------------------------------------
ppi_gene_filtered_df <- ppi_pair_df %>%
  filter(GENE %in% c("IDH2", "MDM4", "PBRM1", "PRKCI", "MECOM", "TET2", "GOLPH3", "FGFR4", "SQSTM1", "RACK1", "ARID1B", "QKI", "OSBPL3", "MYC", "JAK2") | pair_pro %in% c("CECR2:CECR2"))

# annotate degs -----------------------------------------------------------
deg_sample_filtered_df <- merge(deg_sample_df %>%
                                  filter(p_val < 0.05), 
                                ppi_gene_filtered_df, by.x = c("deg_gene_symbol"), by.y = c("SUB_GENE"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, id_aliquot_wu_filter, ".DEG_selected.tsv")
write.table(x = deg_sample_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)


