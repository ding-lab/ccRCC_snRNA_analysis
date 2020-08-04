# Yige Wu @WashU March 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input filtered deg
genes_plot <- c("SLC16A9", "FTCD", "GLDC", "KCNJ15", "PDZD2", "RGN", "ABCC2")
## input median difference
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/compare_bulk_protein_tumor_vs_normal/20200728.v1/Bulk_Protein_Tumor_vs_Normal.Wilcox.20200728.v1.tsv")

# make data matrix for the heatmap ----------------------------------------
plot_data_df <- deg_df %>%
  filter(gene_symbol %in% genes_plot)
rownames(plot_data_df) <- plot_data_df$gene_symbol
plot_data_df <- plot_data_df[genes_plot,]
plot_data_matrix <- as.matrix(plot_data_df[, "meddiff_exp"])
rownames(plot_data_matrix) <- plot_data_df$gene_symbol

# plot --------------------------------------------------------------------
col_fun = colorRamp2(c(-2, 0), c("#a6611a", "white"))
p <- Heatmap(matrix = plot_data_matrix, col = col_fun,
             show_row_names = T,
             cluster_rows = F,
             name = "Bulk Protein\nExpression Difference\n(Tumor-Normal)")
p
## save heatmap
file2write <- paste0(dir_out, "Heatmap", ".pdf")
pdf(file2write, width = 2.8, height = 5)
draw(object = p)
dev.off()
file2write <- paste0(dir_out, "Heatmap", ".png")
png(file2write, width = 500, height = 500, res = 150)
draw(object = p)
dev.off()
