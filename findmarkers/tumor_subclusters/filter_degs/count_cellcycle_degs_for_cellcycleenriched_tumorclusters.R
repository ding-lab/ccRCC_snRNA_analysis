# Yige Wu @WashU Apr 2021

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
## input DEGs
deg_all_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/unite_degs/unite_degs_for_tumor_manualcluster/20210413.v1/TumorManualCluster.DEGs.Wilcox.Minpct0.1.Logfc0.min.diff.pct0.1.tsv")
## input meta data
metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210423.v1/meta_data.20210423.v1.tsv")
## input by cluster enrichment assignment
enrich_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/assign_tumorcluster_by_msigdb_geneset_scores/20210421.v3/MsigDB_Hallmark.Top15GeneSets.4Module.Enrichment.tsv")
## input immune genes
pathway2genes_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/pathway/map_genes2pathway_w_avgexp_sct_data/20210414.v1/TumorCluster.DEG2Pathway.tsv")

# process genes for filtering -------------------------------------------------
cellcycle_genes_df <- pathway2genes_df %>%
  filter(GeneSet_Name %in% c("HALLMARK_MITOTIC_SPINDLE", "HALLMARK_MYC_TARGETS_V1", "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_DNA_REPAIR"))
genes_filter <- unique(cellcycle_genes_df$GeneSymbol)
## process gene category
## map category
gene2cat_df <- cellcycle_genes_df %>%
  group_by(GeneSymbol) %>%
  summarise(Category_sim = paste0(GeneSet_Name, collapse = "|"))

# count -------------------------------------------------------------------
## get immune clusters
enrich_df <- enrich_df %>%
  mutate(cluster_name2 = gsub(x = cluster_name, pattern = "\\.", replacement = "-"))
clusternames_filter <- enrich_df$cluster_name2[enrich_df$Immune]
## map easy id
deg_all_df$sample_id <- deg_all_df$easy_id
deg_all_df$easy_id <- mapvalues(x = deg_all_df$sample_id, from = metadata_df$Sample, to = as.vector(metadata_df$Aliquot.snRNA.WU))
## filter
deg_filtered_df <- deg_all_df %>%
  mutate(cluster_name = paste0(easy_id, "_C", (cluster+1))) %>%
  filter(cluster_name %in% clusternames_filter) %>%
  filter(cellcount_group1 >= 50 & cellcount_group2 >= 50) %>%
  # filter(gene %in% genes_filter) %>%
  filter(p_val_adj < 0.05 & avg_logFC > 0)
## count
count_genes_df <- deg_filtered_df %>%
  group_by(gene) %>%
  summarise(count_clusters = n(), mean_avg_logFC = mean(avg_logFC)) %>%
  arrange(desc(count_clusters), desc(mean_avg_logFC))
count_genes_df$Category_sim <- mapvalues(x = count_genes_df$gene, from = gene2cat_df$GeneSymbol , to = as.vector(gene2cat_df$Category_sim))
count_genes_df$Category_sim[count_genes_df$Category_sim == count_genes_df$gene] <- ""
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Count.Cellcycle_DEGs.tsv")
write.table(x = count_genes_df, file = file2write, quote = F, sep = "\t", row.names = F)

