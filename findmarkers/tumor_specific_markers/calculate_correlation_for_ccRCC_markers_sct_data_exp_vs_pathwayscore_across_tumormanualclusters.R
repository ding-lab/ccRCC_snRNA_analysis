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
## input the average expression calculated (SCT)
avgexp_df <- fread(input = "./Resources/Analysis_Results/average_expression/avgeexp_tumorcells_sct_data_by_manualcluster_rm_doublets_on_katmai/20210413.v1/AverageExpression_ByManualTumorSubcluster.20210413.v1.tsv", data.table = F)
## input cell number per cluster
cellnumber_percluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/count_cellnumber_per_manual_cluster_rm_doublet/20210413.v1/CellNumberPerTumorManualCluster.20210413.v1.tsv")
markers_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/overlap_tumor_vs_pt_DEGs_w_tumor_vs_other_DEGs/20210712.v1/ccRCC_markers.Surface.20210712.v1.tsv")
## input pathway scores
pathscores_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/calculate_msigdb_top_geneset_scores/20210419.v1/MSigDB.Hallmark.tsv")

# preprocess -------------------------------------------------
## identify clusters with sufficient cell number
cluster_pass_df <- cellnumber_percluster_df %>%
  filter(Freq >= 50)%>%
  mutate(colname_exp = gsub(x = id_cluster_uniq,pattern = "\\-", replacement = "."))

## add name for the marker groups
genes_process <- markers_df$Gene
genes_process <- genes_process[!(genes_process %in% c("PIK3CB", "ARHGEF28", "PTGER3", "PARD3", "GNG12", "EFNA5", "SPIRE1", "LIFR", "PKP4", "SORBS1", "PTPRM", "FBXO16", "PAM"))]

# format expression data --------------------------------------------------
exp_filtered_long_df <- avgexp_df %>%
  filter(V1 %in% genes2filter) %>%
  melt() %>%
  mutate(id_bycluster_byaliquot = gsub(x = variable, pattern = "SCT.", replacement = "")) %>%
  dplyr::filter(id_bycluster_byaliquot %in% cluster_pass_df$colname_exp) %>%
  mutate(easyid_column = str_split_fixed(string = id_bycluster_byaliquot, pattern = "_", n = 2)[,1]) %>%
  mutate(cluster_name = str_split_fixed(string = id_bycluster_byaliquot, pattern = "_", n = 2)[,2])

# calculate correlation ---------------------------------------------------
p.value_vec <- NULL
rho_vec <- NULL
gene_vec <- NULL
scorename_vec <- NULL
# scorename_tmp <- "EPITHELIAL_MESENCHYMAL_TRANSITION_Score"
# gene_tmp <- "ABCC3"
for (scorename_tmp in colnames(pathscores_df)[1:20]) {
  for (gene_tmp in genes_process) {
    testdata_df <- merge(x = pathscores_df[, c(scorename_tmp, "cluster_name")],
                         y = exp_filtered_long_df %>%
                           filter(V1 == gene_tmp),
                         by.x = c("cluster_name"), by.y = c("id_bycluster_byaliquot"))
    
    spearman_result <- cor.test(x = testdata_df[, "value"], y = testdata_df[, scorename_tmp], method = "spearman")
    p.value_vec <- c(p.value_vec, spearman_result$p.value)
    rho_vec <- c(rho_vec, spearman_result$estimate)
    gene_vec <- c(gene_vec, gene_tmp)
    scorename_vec <- c(scorename_vec, scorename_tmp)
  }
}
spearman_result_df <- data.frame(scorename = scorename_vec, gene = gene_vec, p.value = p.value_vec, rho = rho_vec)
spearman_result_df$fdr <- p.adjust(p = spearman_result_df$p.value, method = "fdr")
spearman_result_df <- spearman_result_df %>%
  arrange(fdr)
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "ccRCC_marker_expression_vs_pathway_score_bytumorcluster.correlation", ".tsv")
write.table(x = spearman_result_df, file = file2write, sep = "\t", row.names = F, quote = F)
