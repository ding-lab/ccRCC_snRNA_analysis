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
version_tmp <- 3
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input degs
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/by_PT_clusters/findmarker_wilcox_PT_C034911_vs_others_on_katmai/20210812.v1/PT.Byclustergroup.wilcox.logfc.threshold0.min.pct0.1.min.diff.pct0.1.AssayRNA.tsv")
## input EMT related genes
gene2highlight_pos_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/tumor_specific_markers/20210701.v1/ccRCC_cells_specific_DEG_with_surface_annotations_from_3DB.txt")
gene2highlight_neg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_by_celltype/run_bycelltype_bysample/Proximal_tubule/Proximal_tubule_specific_DEG_with_surface_annotations_from_3DB.txt")

# filter the emt genes ----------------------------------------------------
gene2highlight_pos_df <- gene2highlight_pos_df %>%
  arrange(desc(avg_log2FC)) %>%
  head(100)
gene2highlight_neg_df <- gene2highlight_neg_df %>%
  arrange(desc(avg_log2FC)) %>%
  head(100)
genes_hl_pos <- c(gene2highlight_df$Gene, "ITGB8", "PIGR")
genes_hl_pos <- unique(genes_hl_pos)
genes_hl_neg <- gene2highlight_neg_df$Gene

# set plotting parameters -------------------------------------------------
## set y bottom threshold
y_bottom <- -log10(0.05)
## set x limits to distinguish colors
# x_pos <- 1
# x_neg <- -1
x_pos <- 0
x_neg <- 0
## colors
color_right_deep <- RColorBrewer::brewer.pal(n = 12, name = "Paired")[12]
color_right_pale <- RColorBrewer::brewer.pal(n = 7, name = "Set2")[7]
color_left_deep <- RColorBrewer::brewer.pal(n = 6, name = "Dark2")[4]
color_left_pale <- RColorBrewer::brewer.pal(n = 6, name = "Set2")[4]

# make data for plotting --------------------------------------------------
plot_data_df <- deg_df %>%
  dplyr::rename(genesymbol = genesymbol_deg) %>%
  mutate(Log10p_val_adj = -log10(x = p_val_adj)) %>%
  mutate(x_plot = avg_log2FC)
## cap y axis
y_cap <- max(plot_data_df$Log10p_val_adj[!is.infinite(plot_data_df$Log10p_val_adj)])
plot_data_df <- plot_data_df %>%
  mutate(y_plot = ifelse(Log10p_val_adj >= y_cap, y_cap, Log10p_val_adj)) %>%
  mutate(text_gene = ifelse((y_plot >= y_bottom) & ((genesymbol %in% genes_hl_pos) & (x_plot >= x_pos)) | ((genesymbol %in% genes_hl_neg) & (x_plot <= x_neg)), genesymbol, NA))

# plot all markers--------------------------------------------------------------------
## plot
p <- ggplot()
p <- p + geom_vline(xintercept = x_pos, linetype = 2, color = "grey70")
p <- p + geom_vline(xintercept = x_neg, linetype = 2, color = "grey70")
p <- p + geom_point(data = subset(plot_data_df, y_plot < y_bottom), mapping = aes(x = x_plot, y = y_plot), alpha = 0.5, size = 0.5, color = "grey70")
p <- p + geom_point(data = subset(plot_data_df, x_plot < x_pos & x_plot > 0 & y_plot >= y_bottom), mapping = aes(x = x_plot, y = y_plot), alpha = 0.5, size = 0.5, color = "grey70")
p <- p + geom_point(data = subset(plot_data_df, x_plot >= x_pos & y_plot >= y_bottom), mapping = aes(x = x_plot, y = y_plot), alpha = 0.5, size = 0.5, color = color_right_deep)
p <- p + geom_point(data = subset(plot_data_df, x_plot > x_neg & x_plot < 0 & y_plot >= y_bottom), mapping = aes(x = x_plot, y = y_plot), alpha = 0.5, size = 0.5, color = "grey70")
p <- p + geom_point(data = subset(plot_data_df, x_plot <= x_neg & y_plot >= y_bottom), mapping = aes(x = x_plot, y = y_plot), alpha = 0.5, size = 0.5, color = color_left_deep)
# p <- p + scale_color_manual(values = c("FDR<0.05 (up)" = "red", "FDR<0.05 (down)" = "blue", "FDR<0.05" = "black", "FDR>=0.05" = "grey80"))
# p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_gene)),
#                          mapping = aes(x = x_plot, y = y_plot, label = text_gene), color = "black", force = 4, fontface = "bold", segment.alpha = 0.5)
p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_gene) & x_plot > 0),
                         mapping = aes(x = x_plot, y = y_plot, label = text_gene), 
                         color = "black", force = 4, fontface = "italic", segment.alpha = 0.5, size = 5, max.overlaps = Inf, xlim = c(0, 6))
p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_gene) & x_plot < 0),
                         mapping = aes(x = x_plot, y = y_plot, label = text_gene), 
                         color = "black", force = 4, fontface = "italic", segment.alpha = 0.5, segment.size = 0.2,
                         size = 5, max.overlaps = Inf, xlim = c(-6, 0))

p <- p + theme_classic()
p <- p + xlab("Log2(Fold-Change)")
p <- p + ylab("-Log10(P-value-adjusted)")
# p <- p + xlim(c(-4, 4))
p <- p + xlim(c(-6, 6))

p <- p + theme(axis.text = element_text(size = 14, color = "black"),
               axis.title = element_text(size = 15))
p
file2write <- paste0(dir_out, "volcano.", "png")
# png(file2write, width = 800, height = 600, res = 150)
png(file2write, width = 1000, height = 800, res = 150)
print(p)
dev.off()

file2write <- paste0(dir_out, "volcano.", "pdf")
pdf(file2write, width = 8, height = 6, useDingbats = F)
# pdf(file2write, width = 6, height = 5, useDingbats = F)
print(p)
dev.off()

