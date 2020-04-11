# Yige Wu @WashU Apr 2020

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

# input denpendencies -----------------------------------------------------
## input average expression by cell type by aliquot
avgexp_bycelltype_byaliquot_df <- fread(input = "./Resources/Analysis_Results/average_expression/averageexpression_bycelltypeshorter_byaliquot_on_katmai/20200411.v1/averageexpression_bycelltypeshorter.30_aliquot_integration.20200411.v1.tsv", data.table = F)

# make adjustment ---------------------------------------------------------
## get unique filtered genes
genes_filtered <- unique(celltypemarker_filtered_df$gene)
length(genes_filtered)
## filter average expression by filtered genes
avgexp_bycelltype_byaliquot_filtered_df <- avgexp_bycelltype_byaliquot_df %>%
  filter(V1 %in% genes_filtered)
## adjust average expression by 

