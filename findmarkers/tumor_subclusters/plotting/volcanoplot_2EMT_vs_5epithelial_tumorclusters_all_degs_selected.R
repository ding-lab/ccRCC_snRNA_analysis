# Yige Wu @WashU Aug 2021

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
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/findmarkers_selected_EMTclusters_vs_epithelialclusters_katmai/20210924.v3/Selected_2EMTclusters_vs_5Epithelialclusters.logfc.threshold0.min.pct0.1.min.diff.pct0.AssaySCT.tsv")
genes_mesenchymal <- c("SERPINE1", "TGFBI", "VIM", "FN1", "WNT5B", "ITGA5", "JUN", "TWIST1")
# genes_mesenchymal <- c("WNT5B", genes_mesenchymal)
genes_epithelal <- c("LRP2", "ABI3BP", "PTGER3", "FRMD3", "SLC28A1", "SLC6A3", "EPB41LA4", "NFIB", "NFIA", "HNF4A", "HNF4G", "CIT")
genes_filter <- c(genes_mesenchymal, genes_epithelal)
genes_filter

# set plotting parameters -------------------------------------------------
## set y bottom threshold
y_bottom <- -log10(0.05)
## colors
color_right_deep <- RColorBrewer::brewer.pal(n = 12, name = "Set1")[1]
color_left_deep <- RColorBrewer::brewer.pal(n = 6, name = "Set1")[2]

# make data for plotting --------------------------------------------------
plot_data_df <- deg_df %>%
  dplyr::rename(genesymbol = genesymbol_deg) %>%
  mutate(Log10p_val_adj = -log10(x = p_val_adj)) %>%
  mutate(x_plot = ifelse(avg_log2FC < -3, -3,
                         ifelse(avg_log2FC > 3, 3, avg_log2FC)))
## cap y axis
y_cap <- max(plot_data_df$Log10p_val_adj[!is.infinite(plot_data_df$Log10p_val_adj)])
## set x limits to distinguish colors
x_pos <- log2(1.2)
x_neg <- -log2(1.2)
x_pos <- quantile(x = plot_data_df$avg_log2FC, 0.925)
x_neg <- quantile(x = plot_data_df$avg_log2FC, 0.025)
plot_data_df <- plot_data_df %>%
  mutate(y_plot = ifelse(Log10p_val_adj >= y_cap, y_cap, Log10p_val_adj)) %>%
  
  mutate(text_gene = ifelse((y_plot >= y_bottom) & ((genesymbol %in% genes_mesenchymal) & (x_plot >= x_pos)) | ((genesymbol %in% genes_epithelal) & (x_plot <= x_neg)), genesymbol, NA))

# plot all markers--------------------------------------------------------------------
## plot
p <- ggplot()
p <- p + geom_point(data = plot_data_df, mapping = aes(x = x_plot, y = y_plot), alpha = 0.5, size = 0.5)
p <- p + geom_point(data = subset(plot_data_df, y_plot < y_bottom), mapping = aes(x = x_plot, y = y_plot), alpha = 0.5, size = 0.5, color = "grey70")
p <- p + geom_point(data = subset(plot_data_df, x_plot >= 0 & y_plot >= y_bottom), mapping = aes(x = x_plot, y = y_plot), alpha = 0.5, size = 0.5, color = color_right_deep)
p <- p + geom_point(data = subset(plot_data_df, x_plot <= 0 & y_plot >= y_bottom), mapping = aes(x = x_plot, y = y_plot), alpha = 0.5, size = 0.5, color = color_left_deep)
p <- p + scale_color_manual(values = c("FDR<0.05 (up)" = "red", "FDR<0.05 (down)" = "blue", "FDR<0.05" = "black", "FDR>=0.05" = "grey80"))
# p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_gene)),
#                          mapping = aes(x = x_plot, y = y_plot, label = text_gene), color = "black", force = 4, fontface = "italic", segment.alpha = 0.5)
p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_gene) & x_plot > 0),
                         mapping = aes(x = x_plot, y = y_plot, label = text_gene),
                         color = "black", force = 5, fontface = "italic", segment.alpha = 0.5, segment.size = 0.2,
                         size = 5, max.overlaps = Inf, xlim = c(0.5, 4.5))
p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_gene) & x_plot < 0),
                         mapping = aes(x = x_plot, y = y_plot, label = text_gene),
                         color = "black", force = 5, fontface = "italic", segment.alpha = 0.5, segment.size = 0.2,
                         size = 5, max.overlaps = Inf, xlim = c(-3.9, 0))
p <- p + theme_classic()
p <- p + ylim(c(0, 350)) + xlim(c(-3.5, 4.5))
p <- p + xlab("Log2(Fold-Change)")
p <- p + ylab("-Log10(P-value-adjusted)")
p <- p + theme(axis.text = element_text(size = 14, color = "black"),
               axis.title = element_text(size = 15))
p
file2write <- paste0(dir_out, "volcano.", "png")
png(file2write, width = 800, height = 600, res = 150)
print(p)
dev.off()

file2write <- paste0(dir_out, "volcano.", "pdf")
pdf(file2write, width = 6, height = 5, useDingbats = F)
print(p)
dev.off()

