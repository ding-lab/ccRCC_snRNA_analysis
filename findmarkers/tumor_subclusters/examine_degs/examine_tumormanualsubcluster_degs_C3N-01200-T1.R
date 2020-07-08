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
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/findallmarkers_wilcox_tumorcells_by_manualsubcluster_with_cnv_diff_using_individualobj_on_katmai/20200630.v7/Tumormanualsubcluster.withCNVDiff.FindMarkers.Wilcox.Minpct0.1.Logfc0.1.tsv")
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200505.v1/meta_data.20200505.v1.tsv")
## input ppi table
ppi_pair_df <- fread(data.table = F, input = "./Resources/Knowledge/Databases/Protein_Protein_Interactions/protein_pair_table_v2.txt")
## input HIF downstream to look up 
ccrcc_downstream_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/write_ccrcc_genetic_event_downstream_genes/20200302.v2/ccRCC_Genetic_Event_Downstream_Genes.20200302.v2.tsv")
## set sample id
id_aliquot_wu_filter <- "C3N-01200-T1"

# filter by sample --------------------------------------------------------
deg_sample_df <- deg_df %>%
  filter(id_aliquot_wu == id_aliquot_wu_filter) %>%
  unique()
deg_sample_sig_df <- deg_sample_df %>%
  filter(p_val_adj < 0.05) %>%
  filter(ident.1 == "C3N-01200-T1_C3" | ident.2 == "C3N-01200-T1_C3")
deg_sample_sig_df1 <- deg_sample_sig_df %>%
  filter(ident.1 == "C3N-01200-T1_C1" | ident.2 == "C3N-01200-T1_C1")
deg_sample_sig_df2 <- deg_sample_sig_df %>%
  filter(ident.1 == "C3N-01200-T1_C2" | ident.2 == "C3N-01200-T1_C2")
degs_overlap <- intersect(deg_sample_sig_df1$deg_gene_symbol, deg_sample_sig_df2$deg_gene_symbol)
deg_sample_sig_df <- deg_sample_sig_df %>%
  filter(deg_gene_symbol %in% degs_overlap)

# filter ppit -------------------------------------------------------------
ppi_gene_filtered_df <- ppi_pair_df %>%
  filter(GENE %in% c("MDM4", "SETD2", "VHL", "PRKCI", "MECOM", "GOLPH3", "FGFR4", "SQSTM1", "RACK1", "MYC", "JAK2", "HIF1A", "EPAS1"))

# annotate degs -----------------------------------------------------------
deg_sample_filtered_df1 <- merge(deg_sample_sig_df, ppi_gene_filtered_df, 
                                by.x = c("deg_gene_symbol"), by.y = c("SUB_GENE"))
deg_sample_filtered_df2 <- merge(deg_sample_sig_df, 
                                 ccrcc_downstream_df %>%
                                   filter(source_genesymbol %in% c("HIF1A", "EPAS1", "MYC")), 
                                 by.x = c("deg_gene_symbol"), by.y = c("target_genesymbol"))
deg_sample_filtered_df <- rbind(deg_sample_filtered_df1 %>%
                                  select(-pair_pro), 
                                deg_sample_filtered_df2 %>%
                                  rename(GENE = source_genesymbol) %>%
                                  mutate(SUB_GENE.is_TF_downstream = T, SUB_GENE.is_TF_upstream = F, SUB_GENE.is_kinase_substrate = F, SUB_GENE.is_phosphatase_substrate = F, SUB_GENE.is_upstream_kinase = F, SUB_GENE.is_upstream_phosphatase = F, SUB_GENE.is_complex_partner = F))
deg_sample_filtered_df <- unique(deg_sample_filtered_df)
deg_sample_filtered_df  %>% select(deg_gene_symbol) %>% unique()
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, id_aliquot_wu_filter, ".DEG_selected.tsv")
write.table(x = deg_sample_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)

