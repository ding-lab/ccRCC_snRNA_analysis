# Yige Wu @WashU June 2022

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "data.table",
  "stringr",
  "ggpubr",
  "ggplot2"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    install.packages(pkg_name_tmp, dependencies = T)
  }
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg_name_tmp)
  }
  library(package = pkg_name_tmp, character.only = T)
}

# input data --------------------------------------------------------------
## input cell type fraction
scores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/calculate_msigdb_geneset_scores/20210805.v1/MSigDB.Hallmark.tsv")
## input gene sets to test
genesets_test_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/pathway/count_ora_sig_genesets_in_up_degs_across_samples/20220606.v1/Count_gene_set_in_up_tumorcluster_degs.20220606.v1.tsv")
## load CNV fraction in tumor cells
cnv_3state_count_aliquots <- fread("./Resources/Analysis_Results/copy_number/summarize_cnv_fraction/cnv_fraction_in_tumorcells_per_manualcluster/20201207.v1/fraction_of_tumorcells_with_cnv_by_gene_by_3state.per_manualsubcluster.20201207.v1.tsv", data.table = F)
## input known CNV genes
knowncnvgenes_df <- readxl::read_xlsx(path = "./Resources/Knowledge/Known_Genetic_Alterations/Known_CNV.20200528.v1.xlsx", sheet = "Genes")

# pre-process -------------------------------------------------------------
genesets_test <- genesets_test_df$Description
scoregroups_test <- paste0(gsub(x = genesets_test_df$Description, pattern = "HALLMARK_", replacement = ""), "_Score")

# test --------------------------------------------------------------------
gene_tmp <- "SQSTM1"
result_df <- NULL
for (scoregroup_tmp in scoregroups_test) {
  scores_tmp_df <- scores_df[, c("cluster_name", scoregroup_tmp)]
  colnames(scores_tmp_df) <- c("cluster_name", "score.bycluster")
  scores_tmp_df <- scores_tmp_df %>%
    mutate(cluster_name.formatted = gsub(x = cluster_name, pattern = "\\.", replacement = "-"))
  
  ## initiate
  pearson_pvalue_vec <- NULL
  pearson_coef_vec <- NULL
  spearman_pvalue_vec <- NULL
  spearman_rho_vec <- NULL
  for (gene_tmp in knowncnvgenes_df$Gene_Symbol) {
    ## make plot data
    expected_cnv_direction <- knowncnvgenes_df$CNV_Type[knowncnvgenes_df$Gene_Symbol == gene_tmp]
    cnv_frac_tmp_df <- cnv_3state_count_aliquots %>%
      filter(gene_symbol == gene_tmp) %>%
      filter(cna_3state == expected_cnv_direction)
    test_data_comp_df <- merge(x = scores_tmp_df,
                               y = cnv_frac_tmp_df,
                               by.x = c("cluster_name.formatted"), by.y = c("tumor_subcluster"), all.x = T)
    test_data_nonna_df <- test_data_comp_df %>%
      filter(!is.na(Fraction) & !is.na(score.bycluster))
    if (nrow(test_data_nonna_df) >= 5) {
      pearson_result <- cor.test(x = test_data_comp_df$score.bycluster, y = test_data_comp_df$Fraction, method = "pearson")
      pearson_pvalue_vec <- c(pearson_pvalue_vec, pearson_result$p.value)
      pearson_coef_vec <- c(pearson_coef_vec, pearson_result$estimate)
      
      spearman_result <- cor.test(x = test_data_comp_df$score.bycluster, y = test_data_comp_df$Fraction, method = "spearman")
      spearman_pvalue_vec <- c(spearman_pvalue_vec, spearman_result$p.value)
      spearman_rho_vec <- c(spearman_rho_vec, spearman_result$estimate)
    } else {
      pearson_pvalue_vec <- c(pearson_pvalue_vec, NA)
      pearson_coef_vec <- c(pearson_coef_vec, NA)
      spearman_pvalue_vec <- c(spearman_pvalue_vec, NA)
      spearman_rho_vec <- c(spearman_rho_vec, NA)
    }

  }
  result_tmp_df <- data.frame(gene_symbol = knowncnvgenes_df$Gene_Symbol, cnv_direction = knowncnvgenes_df$CNV_Type,
                              pearson_pvalue = pearson_pvalue_vec, pearson_coef = pearson_coef_vec,
                              spearman_pvalue = spearman_pvalue_vec, spearman_rho = spearman_rho_vec)
  result_tmp_df$geneset <- gsub(x = scoregroup_tmp, pattern = "_Score", replacement = "")
  result_df <- rbind(result_df, result_tmp_df)
}
result_df$pearson_fdr <- p.adjust(p = result_df$pearson_pvalue, method = "fdr")
result_df$spearman_fdr <- p.adjust(p = result_df$spearman_pvalue, method = "fdr")
result_df <- result_df %>%
  arrange(pearson_fdr)

# write output ------------------------------------------------------------
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
file2write <- paste0(dir_out, "Tumorcluster_score_vs_CNV_fraction.", run_id, '.tsv')
write.table(x = result_df, file = file2write, quote = F, row.names = F, sep = "\t")


