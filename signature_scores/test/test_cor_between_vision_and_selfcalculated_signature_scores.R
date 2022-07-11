
# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}

# input dependencies ------------------------------------------------------
## input median(?) signature scores per cluster
visionscores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/signature_scores/process_signature_scores/getmedian_visionscore_bysiggeneset_byintrapatientcluster/20220606.v1/median_vision_score.sig_genesets.byintrapatient_tumor_cluster.20220606.v1.tsv")
## input pathway scores
selfscores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/calculate_msigdb_geneset_scores/20210805.v1/MSigDB.Hallmark.tsv")
## input the pathways to plot
genesets_test_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/pathway/count_ora_sig_genesets_in_up_degs_across_samples/20220606.v1/Count_gene_set_in_up_tumorcluster_degs.20220606.v1.tsv")
## input clusters to plotf
clusters_selected_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/plotting/heatmap/heatmap_hallmark38_scores_by_manual_tumorcluster_sct_data_scaled_selected/20220706.v1/intrapatient_tumorclusters.selected.20220706.v1.tsv")


# preprocess --------------------------------------------------------------
genesets_test <- genesets_test_df$Description
genesets_test <- genesets_test[genesets_test %in% colnames(visionscores_df)]
clusters_test <- clusters_selected_df$cluster_name
visionscores_df <- visionscores_df %>%
  mutate(cluster_name = gsub(x = intrapatient_cluster_name, pattern = "\\-", replacement = "."))
rownames(visionscores_df) <- visionscores_df$cluster_name
rownames(selfscores_df) <- selfscores_df$cluster_name

# test by gene set --------------------------------------------------------
coefficients_vec <- NULL
pvalues_vec <- NULL
geneset_tmp <- genesets_test[1]
for (geneset_tmp in genesets_test) {
  geneset_tmp2 <- paste0(gsub(x = geneset_tmp, pattern = "HALLMARK_", replacement = ""), "_Score")
  visionscores_tmp <- visionscores_df[clusters_test,geneset_tmp]
  selfscores_tmp <- selfscores_df[clusters_test, geneset_tmp2]
  pearson_result <- cor.test(x = visionscores_tmp, y = selfscores_tmp)
  pvalues_vec <- c(pvalues_vec, pearson_result$p.value)
  coefficients_vec <- c(coefficients_vec, pearson_result$estimate)
}
pearson_result_df <- data.frame(geneset_name = genesets_test, coefficient.pearson = coefficients_vec, pvalue.pearson = pvalues_vec)
pearson_result_df$fdr.pearson <- p.adjust(p = pearson_result_df$pvalue.pearson, method = "fdr")

# input test --------------------------------------------------------------
visiontest_df <- fread(data.table = F, input = "./Resources/Analysis_Results/signature_scores/compare_signature_scores/wilcox_eachcluster_vs_rest_bygeneset_res1_30ccRCC_tumorcellsreclustered/20220411.v2/wilcox.eachcluster_vs_rest.res1.30ccRCCtumorcellreclustered,20220411.v2.tsv")
genesets_sig <- genesets_test[genesets_test %in% visiontest_df$gene_set[visiontest_df$fdr < 0.05]]
genesets_sig
