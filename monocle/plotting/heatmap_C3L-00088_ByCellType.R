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
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/monocle/deg/filter_pseudotime_degs_3samples/20200728.v2/PT_and_TumorCells_ByCellType.Pseudotim_DE_genes.Combine3Samples.Filtered.20200728.v2.tsv")

# plot--------------------------------------------------------------------
## specify genes to plot
genes_plot <- deg_df$gene_symbol
## write output
p <- plot_pseudotime_heatmap(obj_monocle[genes_plot,],
                        cluster_rows = F,
                        cores = 1,
                        show_rownames = T,
                        return_heatmap = T)
file2write <- paste0(dir_out, "Pseudotie_DEG.png")
save_pheatmap_png(x = p, filename = file2write, width = 800, height = 300, res = 150)
file2write <- paste0(dir_out, "Pseudotie_DEG.pdf")
save_pheatmap_pdf(x = p, filename = file2write, width = 5, height = 2)


