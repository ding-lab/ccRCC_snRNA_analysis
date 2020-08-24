# Yige Wu @WashU Aug 2020

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
## input the variable gene list
variable_genes_df <- fread(input = "./Resources/Analysis_Results/findvariablefeatures/findvariablefeatures_cellgroup_stroma/20200811.v1/findvariablefeatures.cellgroup.Stroma.20200811.v1.tsv", data.table = F)
## input genes in druggable pathways
genes_druggable_df <- readxl::read_excel(path = "./Resources/Knowledge/Gene_Lists/Targetable_Genes.20200814.xlsx")

# intersect ---------------------------------------------------------------
genes_overlap <- intersect(variable_genes_df$gene, genes_druggable_df$genesymbol)
genes_overlap
View(genes_druggable_df[genes_druggable_df$genesymbol %in% genes_overlap,])


