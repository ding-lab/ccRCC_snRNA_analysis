# Yige Wu @ WashU 2019 Dec
## filter for potential predictive markers for VEGF signaling in endothelial cells

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")


# set run id --------------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)

# set output directory ----------------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input result for correlation with VEGFA ---------------------------------
cor_stats_sup_tab <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/proteomics/calculate_cor_with_druggable_protein/20191212.v1/Spearman_Cor_Stats_with.VEGFA.PRO.tsv", data.table = F)

# input deg table from integrated samples and get the list of endothelial specific genes ---------------------------------
deg_tab <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/integrate_seurat_objects/20191021.v1/Renal.DEGs.Pos.txt", data.table = F)
deg_tab_uniq <- deg_tab %>%
  arrange(gene, -avg_logFC)
deg_tab_uniq <- deg_tab_uniq[!duplicated(deg_tab_uniq$gene),]


# input cluster 2 cell type table -----------------------------------------
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - Intergrated.AllCluster2Cell_Type.20191211.v1.tsv", data.table = F)

# filter for endothelial specific genes in the bulk correlation results-----------
## get endothelial cell cluster id
ec_clusters <- cluster2celltype_tab$Cluster[cluster2celltype_tab$Enriched_Cell_Type_Abbr == "EC_Podo"]

## get endothelial specific genes
ec_deg_genes <- deg_tab_uniq$gene[deg_tab_uniq$cluster == ec_clusters]

## filter correlation results to ec_deg_genes
cor_stats_filtered <- cor_stats_sup_tab %>%
  filter(geneB %in% ec_deg_genes) %>%
  filter(genepair_cor_fdr < 0.05)


# specify downstream TFs --------------------------------------------------
downstream

# specify downstream pathways ---------------------------------------------



# input TF table and get targets for specified TFs----------------------------------------------------------
tf_tab <- fread(input = "./Ding_Lab/Projects_Current/TP53_shared_data/resources/PPI/TF_interactions.txt", data.table = F)


