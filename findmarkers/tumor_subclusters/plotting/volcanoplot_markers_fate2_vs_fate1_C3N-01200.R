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
deg_df <- fread(data.table = F, input = "./Resources/snRNA_Processed_Data/Differentially_Expressed_Genes/C3N-01200.Tumor_Segments.Merged.DEGs.monocle.branch2.tsv")

# filter the emt genes ----------------------------------------------------
genes_highlight <- c("EFNA5", "WNT5B", "FN1", "PLCB4", "STK39", "SERPINE1", "PLOD2")

# set plotting parameters -------------------------------------------------
## set y bottom threshold
y_bottom <- -log10(0.05)
## set x limits to distinguish colors
x_pos <- 1.5
x_neg <- -1.5
## colors
color_right_deep <- RColorBrewer::brewer.pal(n = 12, name = "Paired")[6]
color_right_pale <- RColorBrewer::brewer.pal(n = 12, name = "Paired")[5]
color_left_deep <- RColorBrewer::brewer.pal(n = 12, name = "Paired")[2]
color_left_pale <- RColorBrewer::brewer.pal(n = 12, name = "Paired")[1]

# make data for plotting --------------------------------------------------
plot_data_df <- deg_df %>%
  mutate(genesymbol = Gene) %>%
  mutate(Log10p_val_adj = -log10(x = p_val_adj)) %>%
  mutate(avg_log2FC = avg_logFC/log(2)) %>%
  mutate(x_plot = avg_log2FC)
## cap y axis
y_cap <- max(plot_data_df$Log10p_val_adj[!is.infinite(plot_data_df$Log10p_val_adj)])
plot_data_df <- plot_data_df %>%
  mutate(y_plot = ifelse(Log10p_val_adj >= y_cap, y_cap, Log10p_val_adj)) %>%
  mutate(text_gene = ifelse((y_plot >= y_bottom) & ( (x_plot >= x_pos)) | ((x_plot <= x_neg) ), genesymbol, NA)) %>%
  mutate(size_text = ifelse(genesymbol %in% genes_highlight, "highlight", "regular"))

# plot all markers--------------------------------------------------------------------
## plot
p <- ggplot()
p <- p + geom_vline(xintercept = x_pos, linetype = 2)
p <- p + geom_vline(xintercept = x_neg, linetype = 2)
p <- p + geom_point(data = subset(plot_data_df, y_plot < y_bottom), mapping = aes(x = x_plot, y = y_plot), alpha = 0.7, size = 0.5, color = "grey70")
p <- p + geom_point(data = subset(plot_data_df, x_plot < x_pos & x_plot > 0 & y_plot >= y_bottom), mapping = aes(x = x_plot, y = y_plot), alpha = 0.7, size = 0.5, color = color_right_pale)
p <- p + geom_point(data = subset(plot_data_df, x_plot >= x_pos & y_plot >= y_bottom), mapping = aes(x = x_plot, y = y_plot), alpha = 0.7, size = 0.5, color = color_right_deep)
p <- p + geom_point(data = subset(plot_data_df, x_plot > x_neg & x_plot < 0 & y_plot >= y_bottom), mapping = aes(x = x_plot, y = y_plot), alpha = 0.7, size = 0.5, color = color_left_pale)
p <- p + geom_point(data = subset(plot_data_df, x_plot <= x_neg & y_plot >= y_bottom), mapping = aes(x = x_plot, y = y_plot), alpha = 0.7, size = 0.5, color = color_left_deep)
# p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_gene)),
#                          mapping = aes(x = x_plot, y = y_plot, label = text_gene), color = "black", force = 4, fontface = "bold", segment.alpha = 0.5)
p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_gene) & x_plot > 0),
                         mapping = aes(x = x_plot, y = y_plot, label = text_gene, size = size_text), 
                         color = "black", force = 3, fontface = "italic", segment.alpha = 0.3, xlim = c(0, NA))
p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_gene) & x_plot < 0),
                         mapping = aes(x = x_plot, y = y_plot, label = text_gene), 
                         color = "black", force = 2, fontface = "italic", segment.alpha = 0.3, xlim = c(NA, 0))
p <- p + scale_size_manual(values = c("highlight" = 6, "regular" = 4))
p <- p + theme_bw()
p <- p + xlim(c(-4, 5))
p <- p + ylim(c(0, 350))
p <- p + xlab("log2(Fold-Change)")
p <- p + ylab("-Log10(P-value-adjusted)")
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank())
p <- p + theme(legend.position = "none")
# p <- p + theme(axis.line = element_line(colour = "black"))
p

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "volcano.", "png")
png(file2write, width = 800, height = 700, res = 150)
print(p)
dev.off()
file2write <- paste0(dir_out, "volcano.", "pdf")
pdf(file2write, width = 6, height = 5, useDingbats = F)
print(p)
dev.off()


