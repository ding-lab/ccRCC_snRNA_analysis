# Yige Wu @WashU Jul 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(monocle)
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input monocle object
obj_monocle <- readRDS(file = "./Resources/snRNA_Processed_Data/Monocle/outputs/C3L-00088/combined_subset_pseudotime_qval_1e-10.rds")
## input deg table
deg_df <- fread(data.table = F, input = "./Resources/snRNA_Processed_Data/Monocle/outputs/C3L-00088/pseudotime_associated_genes/DE_genes_ModelBy_Pseudotime.txt")
## input cell type marker genes
genes2celltype_df <- fread(data.table = F, input = "./Resources/Knowledge/Kidney_Markers/Gene2CellType_Tab.20200406.v1.tsv")

# specify genes to plot ---------------------------------------------------
table(genes2celltype_df$Cell_Type_Group)
genes_proximaltubule <- genes2celltype_df$Gene[genes2celltype_df$Cell_Type1  == "Proximal tuble" | genes2celltype_df$Cell_Type_Group == "Malignant_Nephron_Epithelium"]
genes_plot <- deg_df$gene_short_name[deg_df$pval < 0.1 & deg_df$gene_short_name %in% genes_proximaltubule]
genes_plot

# plot --------------------------------------------------------------------
colanno_pdata <- data.frame(pData(obj_monocle))
colanno_df <- colanno_df %>%
  select(Pseudotime)
plot_pseudotime_heatmap(obj_monocle[genes_plot,],
                             add_annotation_col = colanno_df,
                             # num_clusters = 3,
                             cores = 1,
                             show_rownames = T)

plot_genes_branched_heatmap(obj_monocle[genes_plot,],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
