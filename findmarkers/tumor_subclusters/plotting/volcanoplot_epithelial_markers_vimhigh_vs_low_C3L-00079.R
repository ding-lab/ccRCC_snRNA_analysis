# Yige Wu @WashU Sep 2020

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
## input degs
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/findmarkers_C3L-00079_tumorlike_vs_tumorcells/20200909.v1/findmarkers_wilcox_tumorlikecells_vs_tumorcells.logfcthreshold0.693147180559945.minpct0.1.mindiffpct0.1.tsv")
## input EMT related genes
emtgenes_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/write_emt_genes/20200915.v1/EMT_Genes.20200915.v1.tsv")
## input cell type markers
gene2celltype_df <- fread(data.table = F, input = "./Resources/Knowledge/Kidney_Markers/Gene2CellType_Tab.20200911.v1.tsv")

# filter the emt genes ----------------------------------------------------
genes_filter <- unique(c(emtgenes_df$hgnc_symbol[emtgenes_df$gene_function == "Epithelial"], gene2celltype_df$Gene[gene2celltype_df$Cell_Type1 %in% c("Proximal tubule", "Tumor cells")]) )
genes_filter

# make data for plotting --------------------------------------------------
plot_data_df <- deg_df %>%
  dplyr::rename(genesymbol = row_name) %>%
  mutate(Log10p_val_adj = -log10(x = p_val_adj)) %>%
  mutate(x_plot = avg_logFC) %>%
  mutate(text_gene = ifelse(genesymbol %in% genes_filter, genesymbol, NA))
## cap y axis
y_cap <- max(plot_data_df$Log10p_val_adj[!is.infinite(plot_data_df$Log10p_val_adj)])
plot_data_df <- plot_data_df %>%
  mutate(y_plot = ifelse(Log10p_val_adj >= y_cap, y_cap, Log10p_val_adj))
## set y bottom threshold
y_bottom <- -log10(0.05)
## set x limits to distinguish colors
x_pos <- log(2)
x_neg <- -log(2)

# plot all markers--------------------------------------------------------------------
## plot
p <- ggplot()
# p <- p + geom_hline(yintercept = y_bottom, linetype = 2)
p <- p + geom_point(data = subset(plot_data_df, x_plot >= x_pos), mapping = aes(x = x_plot, y = y_plot), alpha = 0.5, size = 0.5, color = "red")
p <- p + geom_point(data = subset(plot_data_df, x_plot <= x_neg), mapping = aes(x = x_plot, y = y_plot), alpha = 0.5, size = 0.5, color = "blue")
# p <- p + scale_color_manual(values = c("FDR<0.05 (up)" = "red", "FDR<0.05 (down)" = "blue", "FDR<0.05" = "black", "FDR>=0.05" = "grey80"))
p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_gene)),
                         mapping = aes(x = x_plot, y = y_plot, label = text_gene), color = "black", force = 4, fontface = "bold", segment.alpha = 0.5)
# p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_gene) & avg_logFC <= 0),
#                          mapping = aes(x = avg_logFC, y = y_capped, label = text_gene), color = "black", force = 2)
p <- p + theme_bw()
p <- p + ggtitle(label = paste0("C3L-00079 ", "VIM-high transitional cells vs tumor cells"))
# p <- p + ylim(c(0, y_limit))
p <- p + xlim(c(-3, 3))
p
file2write <- paste0(dir_out, "volcano.", "png")
png(file2write, width = 800, height = 600, res = 150)
print(p)
dev.off()


