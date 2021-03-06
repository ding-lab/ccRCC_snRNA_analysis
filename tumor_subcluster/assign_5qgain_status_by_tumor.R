# Yige Wu @WashU Dec 2020
## 3 categories
### partial 3p  loss: require for each CNV, the cluster with the most %cells with CNV sould be over 10% while the cluster with the least %cells with CNVs should be M10%

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
## load CNV fraction in tumor cells
cnv_3state_count_aliquots <- fread("./Resources/Analysis_Results/copy_number/summarize_cnv_fraction/cnv_fraction_in_tumorcells_per_manualcluster/20201207.v1/fraction_of_tumorcells_with_cnv_by_gene_by_3state.per_manualsubcluster.20201207.v1.tsv", data.table = F)
table(cnv_3state_count_aliquots$tumor_subcluster)
## input known CNV genes
knowncnvgenes_df <- readxl::read_xlsx(path = "./Resources/Knowledge/Known_Genetic_Alterations/Known_CNV.20200528.v1.xlsx", sheet = "Genes")
## load meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200427.v1/meta_data.20200427.v1.tsv", data.table = F)
## set genes to filter
genes_filtered <- c("SQSTM1", "RACK1")

# make data frame for assign cnv type for each tumor --------------------------------------------
cnvfraction_df <- cnv_3state_count_aliquots
## add aliquot.wu
cnvfraction_df$aliquot.wu <- mapvalues(x = cnvfraction_df$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
## add expected cna type
cnvfraction_df$gene_expected_state <- mapvalues(x = cnvfraction_df$gene_symbol, from = knowncnvgenes_df$Gene_Symbol, to = as.vector(knowncnvgenes_df$CNV_Type))
## filter genes
cnvfraction_filtered_df <- cnvfraction_df %>%
  filter(gene_symbol %in% genes_filtered) %>%
  filter(!(grepl(x = tumor_subcluster, pattern = "CNA")))
## get 0 values to add
cnvdetected_bystate_df <- cnvfraction_filtered_df %>%
  select(aliquot.wu, gene_symbol, cna_3state, tumor_subcluster) %>%
  table() %>%
  as.data.frame()
cnvdetected_df <- cnvfraction_filtered_df %>%
  select(aliquot.wu, gene_symbol, tumor_subcluster) %>%
  table() %>%
  as.data.frame()
cnvdetected_bystate_df <- merge(x = cnvdetected_bystate_df %>%
                                  rename(Freq.bygene.bystate = Freq), 
                                y = cnvdetected_df %>%
                                  rename(Freq.bygene = Freq), by = c("aliquot.wu", "gene_symbol", "tumor_subcluster"), all.x = T)
cnvdetected_bystate_filtered_df <- cnvdetected_bystate_df %>%
  filter(cna_3state == "Gain") %>%
  filter(Freq.bygene.bystate == 0 & Freq.bygene > 0)
## filter out copy neutral ones
cnvfraction_merged_df <- cnvfraction_filtered_df %>%
  filter(gene_expected_state == cna_3state) %>%
  mutate(id_aliquot_cluster = tumor_subcluster) 
cnvfraction_merged_df <- merge(x = cnvfraction_merged_df,
                                 y = cnvdetected_bystate_filtered_df,
                                 by = c("aliquot.wu", "gene_symbol", "cna_3state", "tumor_subcluster"), all = T)
## add 0 values
cnvfraction_merged_df$Fraction[is.na(cnvfraction_merged_df$Fraction)] <- 0

# examine each aliquot ----------------------------------------------------
## determine samples with minimal 3p loss
cnvfraction_minimal_df <- cnvfraction_merged_df %>%
  group_by(aliquot.wu) %>%
  summarize(Is_minimal_5qgain = ifelse(all(Fraction < 0.1), TRUE, FALSE))
## determine samples with partial 3p loss
cnvfraction_partial_bygene_df <- cnvfraction_merged_df %>%
  group_by(aliquot.wu, gene_symbol) %>%
  summarize(Is_partial_loss = ifelse(any(Fraction < 0.1) & any(Fraction > 0.1), TRUE, FALSE))
cnvfraction_partial_df <- cnvfraction_partial_bygene_df %>%
  group_by(aliquot.wu) %>%
  summarize(Is_partial_5qgain = ifelse(any(Is_partial_loss), TRUE, FALSE))
## finalize sample assignment
cnvtype_bysample_df <- cnvfraction_minimal_df %>%
  select(aliquot.wu) %>%
  mutate(Group_5qgain = ifelse(aliquot.wu %in% cnvfraction_minimal_df$aliquot.wu[cnvfraction_minimal_df$Is_minimal_5qgain], "Minimal",
                               ifelse(aliquot.wu %in% cnvfraction_partial_df$aliquot.wu[cnvfraction_partial_df$Is_partial_5qgain], "Partial", "Complete"))) %>%
  arrange(Group_5qgain)
# > table(cnvtype_bysample_df$Group_5qgain)
# Complete  Minimal  Partial 
# 14        7        8 

# count clusters with/without 3p loss -------------------------------------
## get genes in samples that have minimal 3p loss in all clusters
genes_bysample_minimal_df <- cnvfraction_merged_df %>%
  group_by(gene_symbol, aliquot.wu) %>%
  summarise(Is_gene_minimal_in_sample = all(Fraction < 0.1))
cnvfraction_merged_df <- merge(x = cnvfraction_merged_df,
                               y = genes_bysample_minimal_df,
                               by = c("gene_symbol", "aliquot.wu"), all.x = T)
cnvfraction_merged_df$Group_5qgain <- mapvalues(x = cnvfraction_merged_df$aliquot.wu, from = cnvtype_bysample_df$aliquot.wu, to = as.vector(cnvtype_bysample_df$Group_5qgain))
cnvfraction_subclonal_df <- cnvfraction_merged_df %>%
  filter(!Is_gene_minimal_in_sample) %>%
  filter(Fraction < 0.1) %>%
  filter(Group_5qgain == "Partial") %>%
  select(tumor_subcluster) %>%
  unique()
nrow(cnvfraction_subclonal_df)
cnvfraction_merged_df %>%
  filter(Group_5qgain == "Partial") %>%
  select(tumor_subcluster) %>%
  unique() %>%
  nrow()
cnvfraction_merged_df %>%
  filter(Group_5qgain == "Complete") %>%
  select(tumor_subcluster) %>%
  unique() %>%
  nrow()
cnvfraction_merged_df %>%
  filter(Group_5qgain == "Minimal") %>%
  select(tumor_subcluster) %>%
  unique() %>%
  nrow()

count_subclonal_df <- cnvfraction_merged_df %>%
  filter(!Is_gene_minimal_in_sample) %>%
  filter(Fraction < 0.1) %>%
  select(tumor_subcluster, Group_5qgain) %>%
  unique() %>%
  group_by(Group_5qgain) %>%
  summarise(number_clusters = n())


# write outputs -----------------------------------------------------------
file2write <- paste0(dir_out, "5qGain_Type_Assignment_Per_Tumor.", run_id, ".tsv")
write.table(x = cnvtype_bysample_df, file = file2write, quote = F, row.names = F, sep = "\t")
file2write <- paste0(dir_out, "5qGain_Type_Assignment_Per_TumorCluster.", run_id, ".tsv")
write.table(x = cnvtype_bysample_df, file = file2write, quote = F, row.names = F, sep = "\t")
