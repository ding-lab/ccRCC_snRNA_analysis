# Yige Wu @WashU Feb 2020

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
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_with_cna/findallmarkers_genelevel_expected_cnv_vs_neutral_in_tumorcells/20200518.v1/ExpectedCNV_vs_Neutral..FindAllMarkers.Wilcox..20200518.v1.tsv")
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200505.v1/meta_data.20200505.v1.tsv")
## input ppi table
ppi_pair_df <- fread(data.table = F, input = "./Resources/Knowledge/Databases/Protein_Protein_Interactions/protein_pair_table_v2.txt")
## set cnv gene to examine
gene_cna <- "VHL"
ids2examine <- "C3N-01200-T1"

# add readable id and filter by cnv gene ----------------------------------
## add id
deg_df$id_aliquot_wu <- mapvalues(x = deg_df$id_aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
## filter by cnv gene
deg_filtered_df <- deg_df %>%
  filter(cna_gene_symbol == gene_cna) %>%
  filter(id_aliquot_wu %in% ids2examine) %>%
  arrange(p_val_adj)
deg_filtered_df %>%
  filter(p_val_adj < 0.05) %>%
  nrow()

# examine ubiquitin_proteasome_genes ---------------------------------------------
deg_filtered_df1 <- deg_filtered_df %>%
  filter(de_gene_symbol %in% ubiquitin_proteasome_genes)
# VHL itself is not differentially expressed, may due to the pct expressed threshold during deg analysis
# examine vhl interacting genes ---------------------------------------------
deg_filtered_df2 <- merge(deg_filtered_df, ppi_pair_df, by.x = c("cna_gene_symbol", "de_gene_symbol"), by.y = c("GENE", "SUB_GENE"))
deg_filtered_df2 <- deg_filtered_df2 %>%
  arrange(p_val_adj)
deg_filtered_df2 %>%
  filter(p_val_adj < 0.05) %>%
  nrow()

table(data.frame(deg_filtered_df$p_val_adj < 0.05, deg_filtered_df$de_gene_symbol %in% deg_filtered_df2$de_gene_symbol)) %>% fisher.test()
## VHL interacting genes are not enriched in significant degs

# examine HIF interacting genes ---------------------------------------------
deg_filtered_df3 <- merge(deg_filtered_df %>%
                            mutate(GENE = "HIF1A"), ppi_pair_df, by.x = c("GENE", "de_gene_symbol"), by.y = c("GENE", "SUB_GENE"))
deg_filtered_df3 <- deg_filtered_df3 %>%
  arrange(p_val_adj)

deg_filtered_df4 <- merge(deg_filtered_df %>%
                            mutate(GENE = "EPAS1"), ppi_pair_df, by.x = c("GENE", "de_gene_symbol"), by.y = c("GENE", "SUB_GENE"))
deg_filtered_df4 <- deg_filtered_df4 %>%
  arrange(p_val_adj)

