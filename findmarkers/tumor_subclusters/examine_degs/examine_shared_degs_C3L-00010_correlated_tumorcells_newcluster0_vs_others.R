# Yige Wu @WashU Aug 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(clusterProfiler)
library(biomaRt)
library(org.Hs.eg.db)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input degs shared by all cluster comparisons
deg_shared_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/examine_degs/intersect_degs_C3L-00010_correlated_tumorcells_newcluser0_vs_others/20200818.v1/shared_degs_C3L-00010_correlated_tumorcells_newcluster0vsothers20200818.v1.tsv")
## input genes in oncogenic process (heterogeneity related)
genes_oncogenic_process_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/write_markers_by_intratumorheterogeneity_types/20200504.v1/markergenes_by_intratumorheterogeneity_types.20200504.v1.tsv")
## input genes in druggable pathways
genes_druggable_df <- readxl::read_excel(path = "./Resources/Knowledge/Gene_Lists/Targetable_Genes.20200814.xlsx")

# examine SMGs ----------------------------------------------------------------
View(deg_shared_df[deg_shared_df$gene %in% ccRCC_SMGs,])

# examine genes in pathgenetic pathways ----------------------------------------------------------------
View(deg_shared_df[deg_shared_df$gene %in% genes_pathogenicpathways_in_ccRCC,])

# examine genes in druggable pathways ----------------------------------------------------------------
genes_overlap <- intersect(deg_shared_df$gene, genes_oncogenic_process_df$gene_symbol)
View(deg_shared_df[deg_shared_df$gene %in% genes_overlap,])
View(genes_oncogenic_process_df[genes_oncogenic_process_df$gene_symbol %in% genes_overlap,])

# examine genes in druggable pathways ----------------------------------------------------------------
genes_overlap <- intersect(deg_shared_df$gene, genes_druggable_df$genesymbol)
View(deg_shared_df[deg_shared_df$gene %in% genes_overlap,])
