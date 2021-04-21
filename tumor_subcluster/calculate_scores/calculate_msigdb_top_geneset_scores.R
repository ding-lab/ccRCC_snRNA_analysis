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
## input the average expression calculated (SCT)
avgexp_df <- fread(input = "./Resources/Analysis_Results/average_expression/avgeexp_tumorcells_sct_data_by_manualcluster_rm_doublets_on_katmai/20210413.v1/AverageExpression_ByManualTumorSubcluster.20210413.v1.tsv", data.table = F)
## input the genes to plot
pathway2genes_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/pathway/ora_msigdb_tumor_manualsubcluster_up_degs/20210413.v1/ORA_Results.tsv")
## input meta data
id_metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210322.v1/meta_data.20210322.v1.tsv")
## input cell number per cluster
cellnumber_percluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/count_cellnumber_per_manual_cluster_rm_doublet/20210413.v1/CellNumberPerTumorManualCluster.20210413.v1.tsv")

# preprocess --------------------------------------------------------------
## add id for mapping
pathway2genes_df$sample_id = pathway2genes_df$easy_id
pathway2genes_df$easy_id <- mapvalues(x = pathway2genes_df$sample_id, from = id_metadata_df$Sample, to = as.vector(id_metadata_df$Aliquot.snRNA.WU))
pathway2genes_df <- pathway2genes_df %>%
  mutate(colname_exp = paste0(gsub(x = easy_id,pattern = "\\-", replacement = "."), "_C", (cluster+1)))
## identify clusters with sufficient cell number
cluster_pass_df <- cellnumber_percluster_df %>%
  filter(Freq >= 50)%>%
  mutate(colname_exp = gsub(x = id_cluster_uniq,pattern = "\\-", replacement = "."))

# specify genes to filter -------------------------------------------------
genesetnames_plot <- c("HALLMARK_MITOTIC_SPINDLE", "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_UV_RESPONSE_DN", "HALLMARK_ALLOGRAFT_REJECTION",
                       "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_HYPOXIA", "HALLMARK_MTORC1_SIGNALING",
                       "HALLMARK_COMPLEMENT", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                       "HALLMARK_KRAS_SIGNALING_UP", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_DNA_REPAIR", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_IL2_STAT5_SIGNALING",
                       "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_PI3K_AKT_MTOR_SIGNALING", "HALLMARK_TGF_BETA_SIGNALING", "HALLMARK_ESTROGEN_RESPONSE_EARLY")
## add name for the marker groups
pathway2genes_filtered_df <- pathway2genes_df %>%
  filter(p.adjust < 0.05) %>%
  filter(Description %in% genesetnames_plot)

pathway2genes_list <- sapply(pathway2genes_filtered_df$geneID, function(x) {
  genes_vec <- str_split(string = x, pattern = "\\/")[[1]]
  return(genes_vec)
})
genes2filter <- unique(unlist(pathway2genes_list))
gene2pathway_df <- data.frame(GeneSet_Name = rep(x = pathway2genes_filtered_df$Description, sapply(X = pathway2genes_list, FUN = function(x) {
  length_vec <- length(x)
  return(length_vec)
})), GeneSymbol = unlist(pathway2genes_list))
gene2pathway_df <- gene2pathway_df %>%
  filter(GeneSymbol %in% avgexp_df$V1)

# format expression data --------------------------------------------------
plot_data_long_df <- avgexp_df %>%
  filter(V1 %in% genes2filter) %>%
  melt() %>%
  mutate(id_bycluster_byaliquot = gsub(x = variable, pattern = "SCT.", replacement = "")) %>%
  filter(id_bycluster_byaliquot %in% cluster_pass_df$colname_exp) %>%
  mutate(easyid_column = str_split_fixed(string = id_bycluster_byaliquot, pattern = "_", n = 2)[,1]) %>%
  mutate(cluster_name = str_split_fixed(string = id_bycluster_byaliquot, pattern = "_", n = 2)[,2])
## filter out non-tumor and NA tumor cluster
plot_data_long_df <- plot_data_long_df %>%
  filter(!(cluster_name %in% c("", "CNA")))
## make matrix
plot_data_wide_df <- dcast(data = plot_data_long_df, formula = V1 ~ id_bycluster_byaliquot, value.var = "value")
plot_data_raw_mat <- as.matrix(plot_data_wide_df[,-1])
## add row names
rownames(plot_data_raw_mat) <- plot_data_wide_df$V1
## scale by row
plot_data_mat <- t(apply(plot_data_raw_mat, 1, scale))
rownames(plot_data_mat) <- rownames(plot_data_raw_mat)
colnames(plot_data_mat) <- colnames(plot_data_raw_mat)

# calculate geneset score -------------------------------------------------
colanno_df <- NULL
for (geneset_tmp in genesetnames_plot) {
  genes_tmp <- gene2pathway_df$GeneSymbol[gene2pathway_df$GeneSet_Name == geneset_tmp]
  score_vec <- colMeans(plot_data_mat[genes_tmp,])*100
  anno_name_tmp <- paste0(gsub(x = geneset_tmp, pattern = "HALLMARK_", replacement = ""), "_Score")
  colanno_tmp_df <- data.frame(score_vec); colnames(colanno_tmp_df) <- anno_name_tmp
  if (is.null(colanno_df)) {
    colanno_df <- colanno_tmp_df
  } else {
    colanno_df <- cbind(colanno_df, colanno_tmp_df)
  }
}
colanno_df$cluster_name <- rownames(colanno_df)
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "MSigDB.Hallmark.tsv")
write.table(x = colanno_df, file = file2write, quote = F, sep = "\t", row.names = F)
