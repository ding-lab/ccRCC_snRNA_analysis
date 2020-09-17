# Yige Wu @WashU Sep 2020

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
## input degs to TF table
deg2tf_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/annotate_deg/annotate_tumor_vs_normal_degs_to_tf_by_snatac/20200916.v1/DEGs_with_DA_peaks.20200916.v1.tsv")
## input motif to gene symbol
motif2tfgene_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/annotate_deg/annotate_tumor_vs_normal_degs_to_tf_by_snatac/20200916.v1/Motif_to_GeneSymbol.20200916.v1.tsv")
## specify top motifs
n_top <- 6
## specify the cell type
celltype_process <- "Tumor cells"

# get the genes to filter -------------------------------------------------
genes_tf_filter <- unique(motif2tfgene_df$source_genesymbol[motif2tfgene_df$rank_motif_by_avglogFC <= n_top & motif2tfgene_df$Cell_type.filename == celltype_process])
genes_tf_filter
## choose TF without too much in common
genes_tf_filter <- c("NR3C1", "REL", "NFIC", "HIF1A")

# filter to top 5 motifs --------------------------------------------------
deg2tf_filtered_df <- deg2tf_df %>%
  filter(Cell_type.filename == "Tumor cells") %>%
  filter(source_genesymbol %in% genes_tf_filter) %>%
  mutate(source_genesymbol2 = source_genesymbol)
unique(deg2tf_filtered_df$target_genesymbol)

# write gene annotation ---------------------------------------------------
gene_anno_df <- data.frame(name = unique(c(deg2tf_filtered_df$source_genesymbol, deg2tf_filtered_df$target_genesymbol)))
gene_anno_df <- gene_anno_df %>%
  mutate(gene_type = ifelse(name %in% deg2tf_filtered_df$source_genesymbol, "TF", "TF Target")) %>%
  mutate(gene_color = ifelse(gene_type == "TF", as.character(name), ""))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "tumor_degs_annotated_with_top_motifs.", run_id, ".tsv")
write.table(x = deg2tf_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "genes_annotation.", run_id, ".tsv")
write.table(x = gene_anno_df, file = file2write, quote = F, sep = "\t", row.names = F)

