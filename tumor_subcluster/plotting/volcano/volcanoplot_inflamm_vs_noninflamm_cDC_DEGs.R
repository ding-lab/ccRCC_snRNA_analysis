# Yige Wu @WashU May 2022

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "ggplot2",
  "ggrastr",
  "ggrepel"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    print(paste0(pkg_name_tmp, "is being installed!"))
    BiocManager::install(pkgs = pkg_name_tmp, update = F)
    install.packages(pkg_name_tmp, dependencies = T)
  }
  print(paste0(pkg_name_tmp, " is installed!"))
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input -------------------------------------------------------------------
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_by_celltype/findmarker_wilcox_cDC_inflamm_tumor_vs_noninflamm_katmai/20220615.v1/logfcthreshold.0.minpct.0.mindiffpct.0.tsv")
## input genes to highlight
immune_genes_df1 <- readxl::read_xlsx(path = "./Resources/Knowledge/Gene_Lists/Immune/Cytokines/Cytokine_Genes.20210420.xlsx")
immune_genes_df2 <- fread(data.table = F, input = "./Resources/Knowledge/Gene_Lists/Immune/Cell_state_markers.txt")

# set plotting parameters -------------------------------------------------
color_red <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[1]
color_blue <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[2]
p_val_adj_sig_cutoff <- 0.05
genes_highlight <- unique(c(immune_genes_df1$gene_symbol, immune_genes_df2$Gene))
genes_highlight <- unique(c(immune_genes_df1$gene_symbol))
# genes_highlight <- unique(c(immune_genes_df2$Gene))

# make data for plotting --------------------------------------------------
text_up <- paste0("Up (", length(which(deg_df$p_val_adj < 0.05 & deg_df$avg_log2FC > 0)), ")")
text_down <- paste0("Down (", length(which(deg_df$p_val_adj < 0.05 & deg_df$avg_log2FC < 0)), ")")

plot_data_df <- deg_df %>%
  rename(genesymbol_deg = gene_symbol) %>%
  mutate(log10p_val_adj = -log10(p_val_adj)) %>%
  mutate(foldchange_type = ifelse(p_val_adj < p_val_adj_sig_cutoff, ifelse(avg_log2FC > 0, text_up, text_down), "insignificant"))
summary(plot_data_df$log10p_val_adj)
summary(plot_data_df$avg_log2FC)
table(plot_data_df$foldchange_type)
plot_data_df <- plot_data_df %>%
  # mutate(x_plot = ifelse(avg_log2FC < -3, -3, ifelse(avg_log2FC > 3, 3,  avg_log2FC))) %>%
  mutate(x_plot = avg_log2FC) %>%
  mutate(y_plot = ifelse(log10p_val_adj > 350, 350, log10p_val_adj)) %>%
  arrange(desc(foldchange_type))
summary(plot_data_df$avg_log2FC[plot_data_df$p_val_adj < 0.05 & plot_data_df$avg_log2FC > 0])
x_pos_cutoff <- min(plot_data_df$avg_log2FC[plot_data_df$p_val_adj < 0.05 & plot_data_df$avg_log2FC > 0])
x_neg_cutoff <- max(plot_data_df$avg_log2FC[plot_data_df$p_val_adj < 0.05 & plot_data_df$avg_log2FC < 0])
## make colors
colors_deg <- c(color_red, color_blue, "grey50")
names(colors_deg) <- c(text_up, text_down, "insignificant")
# plot all markers size not scaled--------------------------------------------------------------------
## plot
p <- ggplot()
p <- p + geom_vline(xintercept = x_pos_cutoff, linetype = 2, color = "grey70")
p <- p + geom_vline(xintercept = x_neg_cutoff, linetype = 2, color = "grey70")
p <- p + geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey70")
p <- p + geom_point_rast(data = subset(plot_data_df, foldchange_type == "insignificant"), 
                         mapping = aes(x = x_plot, y = y_plot, color = foldchange_type), alpha = 0.5, shape = 16, size = 0.5)
p <- p + geom_point_rast(data = subset(plot_data_df, foldchange_type != "insignificant"), 
                         mapping = aes(x = x_plot, y = y_plot,  color = foldchange_type), alpha = 0.5, shape = 16, size = 0.5)
p <- p + scale_color_manual(values = colors_deg)
p <- p + geom_text_repel(data = subset(plot_data_df, p_val_adj < 0.05 & avg_log2FC > 0 & (genesymbol_deg %in% genes_highlight)),
                         mapping = aes(x = x_plot, y = y_plot, label = genesymbol_deg),
                         color = "black", alpha = 1, size = 4.5, #fontface = "bold",
                         segment.size = 0.2, segment.alpha = 0.5, #min.segment.length = 0,
                         xlim = c(0, NA), max.overlaps = Inf)
p <- p + geom_text_repel(data = subset(plot_data_df, p_val_adj < 0.05 &  (avg_log2FC < 0) & (genesymbol_deg %in% genes_highlight)),
                         mapping = aes(x = x_plot, y = y_plot, label = genesymbol_deg),
                         color = "black", alpha = 1, size = 4.5, #fontface = "bold",
                         segment.size = 0.2, segment.alpha = 0.5, #min.segment.length = 0,
                         xlim = c(NA, 0), max.overlaps = Inf, force = 5)
p <- p + theme_classic()
# p <- p + xlab("Log2(fold change) of tumor-cell sn gene expression for\nBAP1-mutant tumor vs. other tumors(n=19)")
p <- p + xlab("Log2(fold change)")
p <- p + ylab("-Log10(adjusted P value)")
# p <- p + xlim(c(-3, 3))
p <- p + guides(color = guide_legend(title = paste0("DEG direction"), 
                                     title.position = "top", title.theme = element_text(size = 14), 
                                     nrow = 4, override.aes = aes(size = 3), label.theme = element_text(size = 14)))
p <- p + theme(axis.text = element_text(size = 14, color = "black"),
               axis.title = element_text(size = 14),
               legend.position = "right",
               # legend.position = c(0.80, 0.20), legend.box.background = element_blank(), legend.background = element_blank(),
               legend.box = "horizontal")
p
file2write <- paste0(dir_out, paste0(head(genes_highlight, 20), collapse = "_"), ".pdf")
pdf(file2write, width = 7, height = 4, useDingbats = F)
print(p)
dev.off()

