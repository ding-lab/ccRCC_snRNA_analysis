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
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input monocle object
obj_monocle <- readRDS(file = "./Resources/snRNA_Processed_Data/Monocle/outputs/C3L-00088/combined_subset_pseudotime_qval_1e-10.rds")
## input deg table
deg_df <- fread(data.table = F, input = "./Resources/snRNA_Processed_Data/Monocle/outputs/C3L-00088/pseudotime_associated_genes/branch_determining_genes.txt")
## input other interesting genes
genes_interest_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/write_markers_by_intratumorheterogeneity_types/20200504.v1/markergenes_by_intratumorheterogeneity_types.20200504.v1.tsv")
## input cell type marker genes
genes2celltype_df <- fread(data.table = F, input = "./Resources/Knowledge/Kidney_Markers/Gene2CellType_Tab.20200406.v1.tsv")

# specify genes to plot ---------------------------------------------------
genes_proximaltubule <- genes2celltype_df$Gene[genes2celltype_df$Cell_Type1  == "Proximal tuble" | genes2celltype_df$Cell_Type_Group == "Malignant_Nephron_Epithelium"]

# plot for branch 1--------------------------------------------------------------------
branch_plot <- 1
## specify genes to plot
genes_plot <- deg_df$gene_short_name[deg_df$branch_point == branch_plot & deg_df$qval < 0.1 & deg_df$gene_short_name %in% c(genes_proximaltubule, genes_interest_df$gene_symbol)]
genes_plot <- unique(genes_plot)
## write output
file2write <- paste0(dir_out, "Branch", branch_plot, "_CellTypeMarkerChange.png")
png(filename = file2write, width = 1000, height = 1200, res = 150)
plot_genes_branched_heatmap(obj_monocle[genes_plot,],
                            branch_point = branch_plot,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
dev.off()

# plot for branch 2--------------------------------------------------------------------
branch_plot <- 2
## specify genes to plot
genes_plot <- deg_df$gene_short_name[deg_df$branch_point == branch_plot & deg_df$qval < 0.1 & deg_df$gene_short_name %in% c(genes_proximaltubule, genes_interest_df$gene_symbol)]
genes_plot <- unique(genes_plot)
## write output
file2write <- paste0(dir_out, "Branch", branch_plot, "_CellTypeMarkerChange.png")
png(filename = file2write, width = 1000, height = 800, res = 150)
plot_genes_branched_heatmap(obj_monocle[genes_plot,],
                            branch_point = branch_plot,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
dev.off()
