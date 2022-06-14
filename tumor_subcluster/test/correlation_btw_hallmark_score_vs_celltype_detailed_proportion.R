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
## input cell type proportion
celltype_frac_long_df <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/count_fraction/count_celltypedetailed_fraction_per_sample/20220607.v1/CellGroupBarcodes_Number_and_Fraction_per_Sample20220607.v1.tsv")

# pre-process -------------------------------------------------------------
genesets_test <- genesets_test_df$Description
scoregroups_test <- paste0(gsub(x = genesets_test_df$Description, pattern = "HALLMARK_", replacement = ""), "_Score")
cellgroups_test <- unique(celltype_frac_long_df$Cell_group); cellgroups_test <- cellgroups_test[!(cellgroups_test %in% c("Unknown", "Mixed myeloid/lymphoid"))]; cellgroups_test
celltype_frac_wide_df <- dcast(data = celltype_frac_long_df, formula = Aliquot_WU~Cell_group, value.var = "Frac_CellGroupBarcodes_ByAliquot")
celltype_frac_wide_df[is.na(celltype_frac_wide_df)] <- 0


# test --------------------------------------------------------------------
result_df <- NULL
for (scoregroup_tmp in scoregroups_test) {
  ## make plot data
  scores_tmp_df <- scores_df[, c("cluster_name", scoregroup_tmp)]
  colnames(scores_tmp_df) <- c("cluster_name", "score.bycluster")
  scores_tmp_df <- scores_tmp_df %>%
    mutate(sample_id = str_split_fixed(string = cluster_name, pattern = "_", n = 2)[,1]) %>%
    mutate(Aliquot_WU = gsub(x = sample_id, pattern = "\\.", replacement = "-")) %>%
    group_by(Aliquot_WU) %>%
    summarise(score.bycluster = max(score.bycluster))
  test_data_comp_df <- merge(x = scores_tmp_df,
                        y = celltype_frac_wide_df,
                        by = c("Aliquot_WU"), all.x = T)
  
  ## initiate
  pearson_pvalue_vec <- NULL
  pearson_coef_vec <- NULL
  spearman_pvalue_vec <- NULL
  spearman_rho_vec <- NULL
  for (cellgroup_tmp in cellgroups_test) {
    pearson_result <- cor.test(x = test_data_comp_df$score.bycluster, y = test_data_comp_df[, cellgroup_tmp], method = "pearson")
    pearson_pvalue_vec <- c(pearson_pvalue_vec, pearson_result$p.value)
    pearson_coef_vec <- c(pearson_coef_vec, pearson_result$estimate)
    
    spearman_result <- cor.test(x = test_data_comp_df$score.bycluster, y = test_data_comp_df[, cellgroup_tmp], method = "spearman")
    spearman_pvalue_vec <- c(spearman_pvalue_vec, spearman_result$p.value)
    spearman_rho_vec <- c(spearman_rho_vec, spearman_result$estimate)
  }
  result_tmp_df <- data.frame(cell_group = cellgroups_test, score_group = scoregroup_tmp,
                              pearson_fdr.byscoregroup = p.adjust(p = pearson_pvalue_vec, method = "fdr"),
                              pearson_pvalue = pearson_pvalue_vec, pearson_coef = pearson_coef_vec,
                              spearman_fdr.byscoregroup = p.adjust(p = spearman_pvalue_vec, method = "fdr"),
                              spearman_pvalue = spearman_pvalue_vec, spearman_rho = spearman_rho_vec)
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
file2write <- paste0(dir_out, "Tumorcluster_score_max_per_sample_vs_celltype_detailed_fraction.", run_id, '.tsv')
write.table(x = result_df, file = file2write, quote = F, row.names = F, sep = "\t")


