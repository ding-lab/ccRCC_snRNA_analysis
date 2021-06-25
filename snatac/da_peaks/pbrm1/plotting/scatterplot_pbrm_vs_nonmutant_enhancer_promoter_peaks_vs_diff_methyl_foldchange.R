# Yige Wu @WashU Jun 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_methyl/"
setwd(dir_base)
source("./ccRCC_methyl_analysis/load_pkgs.R")
source("./ccRCC_methyl_analysis/functions.R")
source("./ccRCC_methyl_analysis/variables.R")
source("./ccRCC_methyl_analysis/plotting.R")
library(ggpubr)
library(ggrepel)
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
peaks2probes_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/pbrm1/overlap_degs/overlap_pbrm1_vs_nonmutant_enhancer_promoter_peaks_with_degs/20210625.v1/PBRM1_vs_NonMutant_DAP2DEG.20210625.v1.tsv")

# make plot data ----------------------------------------------------------
plotdata_df <- peaks2probes_df %>%
  filter(!is.na(avg_log2FC.snATAC) & !is.na(avg_log2FC.methyl)) %>%
  filter(!is.na(rho)) %>%
  mutate(log10_methyl_rna_cor_fdr = -log10(fdr))
  # select(avg_log2FC.snATAC, avg_log2FC.methyl, Gene, peak2gene_type, fdr.methyl) %>%
  # unique()

# plot just promoter peaks--------------------------------------------------------------------
p <- ggscatter(data = subset(plotdata_df, peak2gene_type == "Promoter"), x = "avg_log2FC.snATAC", y = "avg_log2FC.methyl", 
               color = "rho", size = "log10_methyl_rna_cor_fdr",
               alpha = 0.6,
               add = "reg.line",  # Add regressin line
               add.params = list(color = "blue", fill = "lightgray", linetype = 2), # Customize reg. line
               conf.int = F # Add confidence interval
)
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey")
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
p <- p + stat_cor(method = "pearson", label.x = 0, label.y = 0)
p <- p + scale_colour_gradient2(low = "blue", high = "red")
# p <- p + scale_color_manual(values = c("Promoter" = brewer.pal(n = 4, name = "Dark2")[2], "Enhancer" = brewer.pal(n = 4, name = "Dark2")[1]))
p <- p + guides(color = guide_legend(nrow = 2, override.aes = aes(size = 3), label.theme = element_text(size = 14)),
                size = guide_legend(nrow = 2, label.theme = element_text(size = 14)))
p <- p + theme(axis.text = element_text(size = 14),
               axis.title = element_text(size = 14),
               legend.position = "bottom", legend.box = "horizontal")
# file2write <- paste0(dir_out, "scatterplot_snATAC_methyl_FC.",".png")
# png(file2write, width = 800, height = 900, res = 150)
# print(p)
# dev.off()
file2write <- paste0(dir_out, "scatterplot_snATAC_methyl_FC.Promoter.",".pdf")
pdf(file2write, width = 5.5, height = 6, useDingbats = F)
print(p)
dev.off()

# plot just promoter peaks, filter methyl probes--------------------------------------------------------------------
p <- ggscatter(data = subset(plotdata_df, peak2gene_type == "Promoter" & fdr < 0.05), x = "avg_log2FC.snATAC", y = "avg_log2FC.methyl", 
               color = "rho", size = "log10_methyl_rna_cor_fdr",
               alpha = 0.6,
               add = "reg.line",  # Add regressin line
               add.params = list(color = "blue", fill = "lightgray", linetype = 2), # Customize reg. line
               conf.int = F # Add confidence interval
)
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey")
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
p <- p + stat_cor(method = "pearson", label.x = 0, label.y = 0)
p <- p + scale_colour_gradient2(low = "blue", high = "red")
# p <- p + scale_color_manual(values = c("Promoter" = brewer.pal(n = 4, name = "Dark2")[2], "Enhancer" = brewer.pal(n = 4, name = "Dark2")[1]))
p <- p + guides(color = guide_legend(nrow = 2, override.aes = aes(size = 3), label.theme = element_text(size = 14)),
                size = guide_legend(nrow = 2, label.theme = element_text(size = 14)))
p <- p + theme(axis.text = element_text(size = 14),
               axis.title = element_text(size = 14),
               legend.position = "bottom", legend.box = "horizontal")
# file2write <- paste0(dir_out, "scatterplot_snATAC_methyl_FC.",".png")
# png(file2write, width = 800, height = 900, res = 150)
# print(p)
# dev.off()
file2write <- paste0(dir_out, "scatterplot_snATAC_methyl_FC.Promoter.ExpCorrelatedProbes.",".pdf")
pdf(file2write, width = 5.5, height = 6, useDingbats = F)
print(p)
dev.off()

# plot just enhancer peaks--------------------------------------------------------------------
p <- ggscatter(data = subset(plotdata_df, peak2gene_type == "Enhancer"), x = "avg_log2FC.snATAC", y = "avg_log2FC.methyl", 
               color = "rho", size = "log10_methyl_rna_cor_fdr",
               alpha = 0.8,
               add = "reg.line",  # Add regressin line
               add.params = list(color = "blue", fill = "lightgray", linetype = 2), # Customize reg. line
               conf.int = F # Add confidence interval
)
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey")
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
p <- p + stat_cor(method = "pearson", label.x = 0, label.y = 0)
p <- p + scale_colour_gradient2(low = "blue", high = "red")
# p <- p + scale_color_manual(values = c("Promoter" = brewer.pal(n = 4, name = "Dark2")[2], "Enhancer" = brewer.pal(n = 4, name = "Dark2")[1]))
# p <- p + guides(color = guide_legend(nrow = 1, override.aes = aes(size = 3), label.theme = element_text(size = 14)))
p <- p + guides(color = guide_legend(nrow = 2, override.aes = aes(size = 3), label.theme = element_text(size = 14)),
                size = guide_legend(nrow = 2, label.theme = element_text(size = 14)))
p <- p + theme(axis.text = element_text(size = 14),
               axis.title = element_text(size = 14),
               legend.position = "bottom", legend.box = "horizontal")
# file2write <- paste0(dir_out, "scatterplot_snATAC_methyl_FC.",".png")
# png(file2write, width = 800, height = 900, res = 150)
# print(p)
# dev.off()
file2write <- paste0(dir_out, "scatterplot_snATAC_methyl_FC.Enhancer.",".pdf")
pdf(file2write, width = 5.5, height = 6, useDingbats = F)
print(p)
dev.off()

# plot just enhancer peaks filter methyl probes--------------------------------------------------------------------
p <- ggscatter(data = subset(plotdata_df, peak2gene_type == "Enhancer" & fdr < 0.05), x = "avg_log2FC.snATAC", y = "avg_log2FC.methyl", 
               color = "rho", size = "log10_methyl_rna_cor_fdr",
               alpha = 0.8,
               add = "reg.line",  # Add regressin line
               add.params = list(color = "blue", fill = "lightgray", linetype = 2), # Customize reg. line
               conf.int = F # Add confidence interval
)
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey")
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
p <- p + stat_cor(method = "pearson", label.x = 0, label.y = 0)
p <- p + scale_colour_gradient2(low = "blue", high = "red")
# p <- p + scale_color_manual(values = c("Promoter" = brewer.pal(n = 4, name = "Dark2")[2], "Enhancer" = brewer.pal(n = 4, name = "Dark2")[1]))
# p <- p + guides(color = guide_legend(nrow = 1, override.aes = aes(size = 3), label.theme = element_text(size = 14)))
p <- p + guides(color = guide_legend(nrow = 2, override.aes = aes(size = 3), label.theme = element_text(size = 14)),
                size = guide_legend(nrow = 2, label.theme = element_text(size = 14)))
p <- p + theme(axis.text = element_text(size = 14),
               axis.title = element_text(size = 14),
               legend.position = "bottom", legend.box = "horizontal")
# file2write <- paste0(dir_out, "scatterplot_snATAC_methyl_FC.",".png")
# png(file2write, width = 800, height = 900, res = 150)
# print(p)
# dev.off()
file2write <- paste0(dir_out, "scatterplot_snATAC_methyl_FC.Enhancer.ExpCorrelatedProbes.",".pdf")
pdf(file2write, width = 5.5, height = 6, useDingbats = F)
print(p)
dev.off()

# plot highlight genes ----------------------------------------------------
plotdata_df <- plotdata_df %>%
  mutate(highlight = abs(avg_log2FC.snATAC) > 1 & abs(avg_log2FC.methyl) > 1)

p <- ggscatter(data = subset(plotdata_df, peak2gene_type == "Promoter"), x = "avg_log2FC.snATAC", y = "avg_log2FC.methyl",
               color = "rho", size = "log10_methyl_rna_cor_fdr",
               add = "reg.line",  # Add regressin line
               add.params = list(color = "grey50", linetype = 2), # Customize reg. line
               conf.int = F # Add confidence interval
)
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey")
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
# p <- p + stat_cor(method = "pearson", label.x = 0, label.y = 0)
p <- p + geom_text_repel(data = subset(x = plotdata_df, highlight == T & peak2gene_type == "Promoter"), 
                         mapping = aes(x = avg_log2FC.snATAC, y = avg_log2FC.methyl, label = Gene), 
                         max.overlaps = Inf, min.segment.length = 0)
p <- p + scale_colour_gradient2(low = "blue", high = "red")
p <- p + guides(color = guide_legend(nrow = 2))
p <- p + theme(axis.text = element_text(size = 14),
               axis.title = element_text(size = 14),
               legend.position = "bottom", legend.box = "horizontal")
# file2write <- paste0(dir_out, "scatterplot_snATAC_methyl_FC." , "1.", "Promoter",".png")
# png(file2write, width = 800, height = 900, res = 150)
# print(p)
# dev.off()
file2write <- paste0(dir_out, "scatterplot_snATAC_methyl_FC.","1.", "Promoter",".pdf")
pdf(file2write, width = 5.5, height = 6, useDingbats = F)
print(p)
dev.off()

# plot highlight genes ----------------------------------------------------
plotdata_df <- plotdata_df %>%
  mutate(highlight = (abs(avg_log2FC.snATAC) > 0.5 & abs(avg_log2FC.methyl) > 1))

p <- ggscatter(data = subset(plotdata_df, peak2gene_type == "Promoter"), x = "avg_log2FC.snATAC", y = "avg_log2FC.methyl",
               color = "rho", size = "log10_methyl_rna_cor_fdr",
               add = "reg.line",  # Add regressin line
               add.params = list(color = "grey50", linetype = 2), # Customize reg. line
               conf.int = F # Add confidence interval
)
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey")
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
# p <- p + stat_cor(method = "pearson", label.x = 0, label.y = 0)
p <- p + geom_text_repel(data = subset(x = plotdata_df, highlight == T & peak2gene_type == "Promoter"), 
                         mapping = aes(x = avg_log2FC.snATAC, y = avg_log2FC.methyl, label = Gene), 
                         max.overlaps = Inf, min.segment.length = 0)
p <- p + scale_colour_gradient2(low = "blue", high = "red")
p <- p + guides(color = guide_legend(nrow = 2))
p <- p + theme(axis.text = element_text(size = 14),
               axis.title = element_text(size = 14),
               legend.position = "bottom", legend.box = "horizontal")
# file2write <- paste0(dir_out, "scatterplot_snATAC_methyl_FC." , "1",".png")
# png(file2write, width = 800, height = 900, res = 150)
# print(p)
# dev.off()
file2write <- paste0(dir_out, "scatterplot_snATAC_methyl_FC.","0.5", "Promoter",".pdf")
pdf(file2write, width = 8, height = 9, useDingbats = F)
# pdf(file2write, width = 5.5, height = 6, useDingbats = F)
print(p)
dev.off()

# plot highlight PTPRJ ----------------------------------------------------
plotdata_df <- plotdata_df %>%
  mutate(highlight = (Gene %in% c("PTPRJ")))

p <- ggscatter(data = subset(plotdata_df, peak2gene_type == "Promoter"), x = "avg_log2FC.snATAC", y = "avg_log2FC.methyl",
               color = "rho", size = "log10_methyl_rna_cor_fdr",
               add = "reg.line",  # Add regressin line
               add.params = list(color = "grey50", linetype = 2), # Customize reg. line
               conf.int = F # Add confidence interval
)
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey")
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
# p <- p + stat_cor(method = "pearson", label.x = 0, label.y = 0)
p <- p + geom_text_repel(data = subset(x = plotdata_df, highlight == T & peak2gene_type == "Promoter"), 
                         mapping = aes(x = avg_log2FC.snATAC, y = avg_log2FC.methyl, label = Gene), 
                         max.overlaps = Inf, min.segment.length = 0)
p <- p + scale_colour_gradient2(low = "blue", high = "red")
p <- p + guides(color = guide_legend(nrow = 2))
p <- p + theme(axis.text = element_text(size = 14),
               axis.title = element_text(size = 14),
               legend.position = "bottom", legend.box = "horizontal")
# file2write <- paste0(dir_out, "scatterplot_snATAC_methyl_FC." , "1",".png")
# png(file2write, width = 800, height = 900, res = 150)
# print(p)
# dev.off()
file2write <- paste0(dir_out, "scatterplot_snATAC_methyl_FC.","PTPRJ",".pdf")
pdf(file2write, width = 5.5, height = 6, useDingbats = F)
print(p)
dev.off()

# plot highlight PTPRJ ----------------------------------------------------
plotdata_df <- plotdata_df %>%
  mutate(highlight = (Gene %in% c("CES3")))

p <- ggscatter(data = subset(plotdata_df, peak2gene_type == "Enhancer"), x = "avg_log2FC.snATAC", y = "avg_log2FC.methyl",
               color = "rho", size = "log10_methyl_rna_cor_fdr",
               add = "reg.line",  # Add regressin line
               add.params = list(color = "grey50", linetype = 2), # Customize reg. line
               conf.int = F # Add confidence interval
)
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey")
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
# p <- p + stat_cor(method = "pearson", label.x = 0, label.y = 0)
p <- p + geom_text_repel(data = subset(x = plotdata_df, highlight == T & peak2gene_type == "Enhancer"), 
                         mapping = aes(x = avg_log2FC.snATAC, y = avg_log2FC.methyl, label = Gene), 
                         max.overlaps = Inf, min.segment.length = 0)
p <- p + scale_colour_gradient2(low = "blue", high = "red")
p <- p + guides(color = guide_legend(nrow = 2))
p <- p + theme(axis.text = element_text(size = 14),
               axis.title = element_text(size = 14),
               legend.position = "bottom", legend.box = "horizontal")
# file2write <- paste0(dir_out, "scatterplot_snATAC_methyl_FC." , "1",".png")
# png(file2write, width = 800, height = 900, res = 150)
# print(p)
# dev.off()
file2write <- paste0(dir_out, "scatterplot_snATAC_methyl_FC.","CES3",".pdf")
pdf(file2write, width = 5.5, height = 6, useDingbats = F)
print(p)
dev.off()


