# Yige Wu @WashU May 2020
## for plotting the fraction of cells with CNV per sample in case per cluster CNV distribution

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
## load meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200427.v1/meta_data.20200427.v1.tsv", data.table = F)
## load CNV fraction in tumor cells
cnv_3state_count_aliquots <- fread("./Resources/Analysis_Results/copy_number/summarize_cnv_fraction/cnv_fraction_in_tumorcells_per_manualcluster/20200622.v1/fraction_of_tumorcells_with_cnv_by_gene_by_3state.per_manualsubcluster.20200622.v1.tsv", data.table = F)
## input known CNV genes
knowncnvgenes_df <- readxl::read_xlsx(path = "./Resources/Knowledge/Known_Genetic_Alterations/Known_CNV.20200528.v1.xlsx", sheet = "Genes")

# add additional info to the CNV data --------------------------------------------
cnvfraction_df <- cnv_3state_count_aliquots
## add aliquot.wu
cnvfraction_df$aliquot.wu <- mapvalues(x = cnvfraction_df$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
## add expected cna type
cnvfraction_df$gene_expected_state <- mapvalues(x = cnvfraction_df$gene_symbol, from = knowncnvgenes_df$Gene_Symbol, to = as.vector(knowncnvgenes_df$CNV_Type))
## filter out copy neutral ones
cnvfraction_df <- cnvfraction_df %>%
  filter(cna_3state != "Neutral") %>%
  filter(gene_expected_state == cna_3state) %>%
  filter(!is.na(tumor_subcluster)) %>%
  mutate(id_aliquot_cluster = paste0(aliquot.wu, "_C", (tumor_subcluster + 1))) 
unique(cnvfraction_df$id_aliquot_cluster)
## filter genes
genes_filtered <- unique(cnvfraction_df$gene_symbol[cnvfraction_df$Fraction > 0.5])
cnvfraction_df <- cnvfraction_df %>%
  filter(gene_symbol %in% genes_filtered)

# get pairs of tumor subclusters within the same sample ----------------------------------------------------
cnvfraction_pair_df <- merge(cnvfraction_df %>%
                               select(-tumor_subcluster), 
                             cnvfraction_df %>%
                               select(-tumor_subcluster), 
                             by = c("aliquot.wu", "aliquot", "gene_symbol", "gene_expected_state"), suffixes = c(".1", ".2"))
## filter
cnvfraction_pair_df <- cnvfraction_pair_df %>%
  filter(id_aliquot_cluster.1 != id_aliquot_cluster.2) %>%
  mutate(fraction_diff = Fraction.1 - Fraction.2)
cnvfraction_pair_filtered_df <- cnvfraction_pair_df %>%
  filter(fraction_diff > 0.5)

# write outputs -----------------------------------------------------------
file2write <- paste0(dir_out, "CNVFractionDifference_between_ManualSubclusterPairs_PerGene_PerTumor.", run_id, ".tsv")
write.table(x = cnvfraction_pair_df, file = file2write, quote = F, row.names = F, sep = "\t")
file2write <- paste0(dir_out, "CNVFractionDifference_between_ManualSubclusterPairs_PerGene_PerTumor.Filtered", run_id, ".tsv")
write.table(x = cnvfraction_pair_filtered_df, file = file2write, quote = F, row.names = F, sep = "\t")

