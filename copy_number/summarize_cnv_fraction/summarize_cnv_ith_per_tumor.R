# Yige Wu @WashU May 2020
## for plotting the fraction of cells with CNV per sample in case per cluster CNV distribution

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
# devtools::install_github("teunbrand/ggh4x")
library(ggh4x)
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
cnv_3state_count_aliquots <- fread("./Resources/Analysis_Results/copy_number/summarize_cnv_fraction/cnv_fraction_in_tumorcells_per_manualcluster/20200505.v1/fraction_of_tumorcells_with_cnv_by_gene_by_3state.per_manualsubcluster.20200505.v1.tsv", data.table = F)
table(cnv_3state_count_aliquots$tumor_subcluster)
## input known CNV genes
knowncnvgenes_df <- readxl::read_xlsx(path = "./Resources/Known_Genetic_Alterations/Known_CNV.20200505.v1.xlsx", sheet = "Genes")

# make data frame for assign cnv type for each tumor --------------------------------------------
cnvfraction_df <- cnv_3state_count_aliquots
## add aliquot.wu
cnvfraction_df$aliquot.wu <- mapvalues(x = cnvfraction_df$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
## add expected cna type
cnvfraction_df$gene_expected_state <- mapvalues(x = cnvfraction_df$gene_symbol, from = knowncnvgenes_df$Gene_Symbol, to = as.vector(knowncnvgenes_df$CNV_Type))
## filter out copy neutral ones
cnvfraction_df <- cnvfraction_df %>%
  filter(cna_3state != "Neutral") %>%
  filter(gene_expected_state == cna_3state) %>%
  mutate(id_aliquot_cluster = paste0(aliquot.wu, "C", tumor_subcluster)) 
## filter genes
genes_filtered <- unique(cnvfraction_df$gene_symbol[cnvfraction_df$Fraction > 0.5])
cnvfraction_df <- cnvfraction_df %>%
  filter(gene_symbol %in% genes_filtered)

# examine each aliquot ----------------------------------------------------
cnvfraction_maxdiff_df <- cnvfraction_df %>%
  group_by(aliquot.wu, gene_symbol, cna_3state) %>%
  summarize(Fraction_maxdiff = (max(Fraction) - min(Fraction)))
cnvfraction_ith_type_df <- cnvfraction_maxdiff_df %>%
  group_by(aliquot.wu) %>%
  summarize(CNV_ITH_Type = ifelse(any(Fraction_maxdiff >= 0.5), "Intra-segment_Heterogeneous", "Intra-segment_Homogenous"))

# write outputs -----------------------------------------------------------
file2write <- paste0(dir_out, "CNV_Fraction_MaxDifferenceAcrossManualSubclusters_PerGene_PerTumor.", run_id, ".tsv")
write.table(x = cnvfraction_maxdiff_df, file = file2write, quote = F, row.names = F, sep = "\t")
file2write <- paste0(dir_out, "CNV_Type_Assignment_Per_Tumor.", run_id, ".tsv")
write.table(x = cnvfraction_ith_type_df, file = file2write, quote = F, row.names = F, sep = "\t")

