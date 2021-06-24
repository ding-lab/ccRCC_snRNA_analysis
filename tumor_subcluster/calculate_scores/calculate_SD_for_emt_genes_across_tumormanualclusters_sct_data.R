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
## input cell number per cluster
cellnumber_percluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/count_cellnumber_per_manual_cluster_rm_doublet/20210413.v1/CellNumberPerTumorManualCluster.20210413.v1.tsv")
## input kidney-specific EMT genes
emt_genes_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/combine_pt_with_emt_markers_all/20200920.v1/Kidney_Specific_EMT_Genes.20200920.v1.tsv")
## input mesenchymal markers
gene2pathway_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/pathway/map_genes2pathway_w_avgexp_sct_data/20210414.v1/TumorCluster.DEG2Pathway.tsv")

# specify genes to filter -------------------------------------------------
## identify clusters with sufficient cell number
cluster_pass_df <- cellnumber_percluster_df %>%
  filter(Freq >= 50)%>%
  mutate(colname_exp = gsub(x = id_cluster_uniq,pattern = "\\-", replacement = "."))

## add name for the marker groups
emt_genes_df <- emt_genes_df %>%
  mutate(Text_Gene_Group = ifelse(Gene_Group2 == "Tumor cells", 
                                  "Tumor-cell\nmarkers", 
                                  ifelse(Gene_Group2 %in% c("Epithelial", "Proximal tubule"),
                                         "Epithelial/\nproximal-tubule\nmarkers", paste0(Gene_Group2, "\nmarkers"))))

genes2filter <- emt_genes_df$Gene

# format expression data --------------------------------------------------
plot_data_long_df <- avgexp_df %>%
  filter(V1 %in% genes2filter) %>%
  melt() %>%
  mutate(id_bycluster_byaliquot = gsub(x = variable, pattern = "SCT.", replacement = "")) %>%
  mutate(easyid_column = str_split_fixed(string = id_bycluster_byaliquot, pattern = "_", n = 2)[,1]) %>%
  mutate(cluster_name = str_split_fixed(string = id_bycluster_byaliquot, pattern = "_", n = 2)[,2])
## filter out non-tumor and NA tumor cluster
plot_data_long_df <- plot_data_long_df %>%
  filter(!(cluster_name %in% c("", "CNA")))
## make matrix
plot_data_wide_df <- dcast(data = plot_data_long_df, formula = V1 ~ id_bycluster_byaliquot, value.var = "value")
plot_data_mat <- as.matrix(plot_data_wide_df[,-1])
rownames(plot_data_mat) <- plot_data_wide_df$V1

# filter genes based on variation -----------------------------------------
sd_bygene_df <- data.frame(SD = apply(plot_data_mat,1, sd, na.rm = TRUE), gene = rownames(plot_data_mat))
sd_bygene_df$Gene_Group2 <- mapvalues(x = sd_bygene_df$gene, from = emt_genes_df$Gene, to = as.vector(emt_genes_df$Gene_Group2))
sd_bygene_df <- sd_bygene_df %>%
  arrange(desc(SD))
genes_plot_mesenchymal <- head(as.vector(sd_bygene_df$gene[sd_bygene_df$Gene_Group2 == "Mesenchymal"]), n = 5)
genes_plot_mesenchymal
# [1] "FN1"   "CDH2"  "VIM"   "FOXC2" "ITGB1"
sd_bygene_epithelial_df <- sd_bygene_df %>%
  filter(Gene_Group2 %in% c("Proximal tubule")) %>%
  arrange(desc(SD))
genes_plot_epithelial <- head(x = as.vector(sd_bygene_epithelial_df$gene), n = 5)
# [1] "GPX3"    "SLC17A3" "SLC13A1" "SLC5A12" "CUBN" 
genes_plot <- c(genes_plot_mesenchymal, genes_plot_epithelial)
plot_data_mat <- plot_data_mat[genes_plot,]

# make mesenchymal/epithelial score score --------------------------------------------------
## make mesenchymal score
genes_mesenchymal_score <- genes_plot_mesenchymal
scores_mesenchymal <- colMeans(plot_data_mat[genes_mesenchymal_score,]); scores_mesenchymal
## make mesenchymal score
genes_epithelial_score <- genes_plot_epithelial
scores_epithelial <- colMeans(plot_data_mat[genes_epithelial_score,]); scores_epithelial
## make threshold
cutoff_mesenchymal <- quantile(x = scores_mesenchymal, probs = 0.9); cutoff_mesenchymal
cutoff_epithelial <- quantile(x = scores_epithelial, probs = 0.1); cutoff_epithelial
## annotate gene table
sd_bygene_df <- sd_bygene_df %>%
  mutate(Used_for_emtscore = ifelse(gene %in% genes_mesenchymal_score, "Mesenchymal_Score",
                                    ifelse(gene %in% genes_epithelial_score, "Epithelial_Score", "")))


# make data frame ---------------------------------------------------------
emtscores_df <- data.frame(id_bycluster_byaliquot = colnames(plot_data_mat),
                           Score_mesenchymal = scores_mesenchymal,
                           Score_epithelial = scores_epithelial)
emtscores_df <- emtscores_df %>%
  mutate(Sample_id = str_split_fixed(string = id_bycluster_byaliquot, pattern = "_", n = 3)[,1]) %>%
  mutate(Cluster_id = str_split_fixed(string = id_bycluster_byaliquot, pattern = "_", n = 3)[,2]) %>%
  select(-id_bycluster_byaliquot) %>%
  select(Sample_id, Cluster_id, Score_mesenchymal, Score_epithelial) %>%
  mutate(EMT_potential = ifelse((scores_mesenchymal >= cutoff_mesenchymal) & (scores_epithelial <= cutoff_epithelial),
                                "High", "Low")) %>%
  arrange(EMT_potential, desc(Score_mesenchymal))

# write outputs -----------------------------------------------------------
file2write <- paste0(dir_out, "EMT_scores_by_manual_tumorcluster.", run_id, ".tsv")
write.table(x = emtscores_df, file = file2write, sep = "\t", row.names = F, quote = F)
file2write <- paste0(dir_out, "EMT_scores_genes.", run_id, ".tsv")
write.table(x = sd_bygene_df, file = file2write, sep = "\t", row.names = F, quote = F)
## save thresholds

