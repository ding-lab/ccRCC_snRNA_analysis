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
genesetnames_plot <- c("HALLMARK_MITOTIC_SPINDLE", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_HYPOXIA", "HALLMARK_MTORC1_SIGNALING",
                       "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_COMPLEMENT", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
                       "HALLMARK_KRAS_SIGNALING_UP", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "HALLMARK_DNA_REPAIR", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_IL2_STAT5_SIGNALING",
                       "HALLMARK_INTERFERON_ALPHA_RESPONSE", "HALLMARK_PI3K_AKT_MTOR_SIGNALING", "HALLMARK_TGF_BETA_SIGNALING")

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

# get gene-pathway map-------------------------------------------------
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
gene2pathway_df <- unique(gene2pathway_df)

# make average normalized expression --------------------------------------------------
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
## get average expression
orig_avgexp_vec <- rowMeans(x = plot_data_raw_mat, na.rm = T)
## map expression
gene2pathway_df$Avgexp.sct.data <- mapvalues(x = gene2pathway_df$GeneSymbol, from = names(orig_avgexp_vec), to = orig_avgexp_vec)
gene2pathway_df$Avgexp.sct.data <- as.numeric(gene2pathway_df$Avgexp.sct.data)

# filter top expressed genes in each pathway ----------------------------------------
gene2pathway_top_df <- gene2pathway_df %>%
  group_by(GeneSet_Name) %>%
  slice_max(order_by = Avgexp.sct.data, n = 10)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "TumorCluster.DEG2Pathway.tsv")
write.table(x = gene2pathway_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "TumorCluster.DEG2Pathway.Top.tsv")
write.table(x = gene2pathway_top_df, file = file2write, quote = F, sep = "\t", row.names = F)

