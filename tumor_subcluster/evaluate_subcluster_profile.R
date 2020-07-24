# Yige Wu @WashU March 2020
## running on local
## for calculating the aliquot-pairwise correlation coefficients for averaged expression of all genes

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
## input id meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv", data.table = F)
## input the CNV fractin per subcluster
cnv_df <- fread("./Resources/Analysis_Results/copy_number/summarize_cnv_fraction/cnv_fraction_in_tumorcells_per_manualcluster/20200622.v1/fraction_of_tumorcells_with_cnv_by_gene_by_3state.per_manualsubcluster.20200622.v1.tsv", data.table = F)
## input known CNV genes
knowncnvgenes_df <- readxl::read_xlsx(path = "./Resources/Knowledge/Known_Genetic_Alterations/Known_CNV.20200528.v1.xlsx", sheet = "Genes")
## input te bulk genomics/methylation events
bulk_sn_omicsprofile_df <- fread(input = "./Resources/Analysis_Results/data_summary/merge_bulk_sn_profiles/20200512.v1/bulk_sn_omics_profile.20200512.v1.tsv", data.table = F)

# process CNV data --------------------------------------------------------
## add aliquot.wu
cnv_df$aliquot.wu <- mapvalues(x = cnv_df$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
cnv_df$case <- mapvalues(x = plot_data_df$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Case))
## add cytoband and expected cna type
cnv_df$gene_cytoband <- mapvalues(x = cnv_df$gene_symbol, from = knowncnvgenes_df$Gene_Symbol, to = as.vector(knowncnvgenes_df$Cytoband))
cnv_df$gene_expected_state <- mapvalues(x = cnv_df$gene_symbol, from = knowncnvgenes_df$Gene_Symbol, to = as.vector(knowncnvgenes_df$CNV_Type))
cnv_df <- cnv_df %>%
  mutate(id_aliquot_cluster = paste0(aliquot.wu, "_C", (tumor_subcluster + 1)))
cnv_plot_df <- cnv_df %>%
  filter(cna_3state == gene_expected_state) %>%
  select(id_aliquot_cluster, gene_symbol, Fraction)
## annotate NA data
cnv_na_df <- cnv_df %>%
  select(id_aliquot_cluster, gene_symbol) %>%
  table() %>%
  as.data.frame() %>%
  filter(Freq == 0) %>%
  mutate(Fraction = 2) %>%
  select(id_aliquot_cluster, gene_symbol, Fraction)
cnv_plot_df <- rbind(cnv_plot_df, cnv_na_df)
## make it wide
cnv_wide_df <- dcast(data = cnv_plot_df, formula = id_aliquot_cluster ~ gene_symbol, value.var = "Fraction", fill = 0)
cnv_wide_df[which(x = cnv_wide_df == 2,arr.ind = T)] <- NA
## add aliquot id
cnv_wide_df <- data.frame(cnv_wide_df)
cnv_wide_df <- cnv_wide_df %>%
  dplyr::mutate(Aliquot_snRNA_WU = str_split_fixed(string = id_aliquot_cluster, pattern = "_", n = 2)[,1])

# merge with bulk omics profile -------------------------------------------
omics_bycluster_df <- merge(bulk_sn_omicsprofile_df, cnv_wide_df, by = c("Aliquot_snRNA_WU"), all.y = T)

# add VHL status ----------------------------------------------------------


# add PBRM1 status ----------------------------------------------------------
# add BAP1 status ----------------------------------------------------------
# add SETD2 status ----------------------------------------------------------


# write output ------------------------------------------------------------


