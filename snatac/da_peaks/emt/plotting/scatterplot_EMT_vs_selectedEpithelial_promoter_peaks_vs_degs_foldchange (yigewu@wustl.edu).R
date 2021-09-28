# Yige Wu @WashU Sep 2021
## BAP1_tumorcells_vs_other_tumorcells

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(ggpubr)
library(ggrepel)
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
peaks2degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/emt/overlap_degs/overlap_EMT_vs_selectedEpithelial_diff_promoter_peaks_with_degs/20210927.v1/EMT_vs_selectedEpithelialClusters.20210927.v1.tsv")

# make plot data ----------------------------------------------------------
plotdata_df <- peaks2degs_df %>%
  filter(!is.na(avg_log2FC.snATAC) & !is.na(avg_log2FC.snRNA)) %>%
  # mutate(avg_log2FC.snATAC = ifelse(is.na(avg_log2FC.snATAC), 0, avg_log2FC.snATAC)) %>%
  # mutate(avg_log2FC.snRNA = ifelse(is.na(avg_log2FC.snRNA), 0, avg_log2FC.snRNA)) %>%
  select(avg_log2FC.snATAC, avg_log2FC.snRNA, Gene, peak) %>%
  unique()
nrow(plotdata_df) ## 1096 gene-peak pairs
length(unique(plotdata_df$Gene)) ## 1030 genes
length(unique(plotdata_df$peak)) ## 1096 promoter peaks

# plot no genes ----------------------------------------------------
p <- ggscatter(data = plotdata_df, x = "avg_log2FC.snATAC", y = "avg_log2FC.snRNA", alpha = 0.8, shape = 16, size = 2.5,
               add = "reg.line", add.params = list(color = "grey50", linetype = 2), # Customize reg. line
               conf.int = F # Add confidence interval
)
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey")
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
p <- p + stat_cor(method = "pearson", label.x = 1, label.y = 1.5, size = 4)
# p <- p + geom_text_repel(data = subset(x = plotdata_df, highlight == T), 
#                          mapping = aes(x = avg_log2FC.snATAC, y = avg_log2FC.snRNA, label = Gene), 
#                          max.overlaps = Inf, size = 6,
#                          segment.size = 0.4, segment.alpha = 0.6, min.segment.length = 0, box.padding = 0.5)
# p <-  p + scale_color_manual(values = colors_peak2genetype)
p <- p + guides(color = guide_legend(nrow = 1, override.aes = aes(size = 3), label.theme = element_text(size = 14)))
p <- p + theme(axis.text = element_text(size = 18, color = "black"),
               axis.title = element_text(size = 14),
               legend.position = "bottom", legend.box = "horizontal")
# file2write <- paste0(dir_out, "scatterplot_snATAC_snRNA_FC." , "1",".png")
# png(file2write, width = 800, height = 900, res = 150)
# print(p)
# dev.off()
file2write <- paste0(dir_out, "scatterplot_snATAC_snRNA_FC", ".pdf")
pdf(file2write, width = 5, height = 5, useDingbats = F)
print(p)
dev.off()

# highlight genes with fold change > 1 ----------------------------------------------------
plotdata_df <- plotdata_df %>%
  mutate(highlight = (avg_log2FC.snATAC >= 1 & avg_log2FC.snRNA >= 1) | (avg_log2FC.snATAC <= -1 & avg_log2FC.snRNA <= -1) | (Gene %in% c("VIM", "FN1", "CDH2", "WNT5B")))
p <- ggscatter(data = plotdata_df, x = "avg_log2FC.snATAC", y = "avg_log2FC.snRNA", alpha = 0.8, shape = 16, size = 2.5,
               add = "reg.line", add.params = list(color = "grey50", linetype = 2), # Customize reg. line
               conf.int = F # Add confidence interval
)
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey")
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
p <- p + stat_cor(method = "pearson", label.x = 1, label.y = 1.5, size = 4)
p <- p + geom_text_repel(data = subset(x = plotdata_df, avg_log2FC.snATAC >= 1 & avg_log2FC.snRNA >= 1),
                         mapping = aes(x = avg_log2FC.snATAC, y = avg_log2FC.snRNA, label = Gene),
                         max.overlaps = Inf, size = 5,
                         segment.size = 0.4, segment.alpha = 0.6, min.segment.length = 0, box.padding = 0.5, color = "red")
p <- p + geom_text_repel(data = subset(x = plotdata_df, avg_log2FC.snATAC <= -1 & avg_log2FC.snRNA <= -1),
                         mapping = aes(x = avg_log2FC.snATAC, y = avg_log2FC.snRNA, label = Gene),
                         max.overlaps = Inf, size = 5,
                         segment.size = 0.4, segment.alpha = 0.6, min.segment.length = 0, box.padding = 0.5, color = "blue", xlim = c(-4, -1))
p <- p + guides(color = guide_legend(nrow = 1, override.aes = aes(size = 3), label.theme = element_text(size = 14)))
p <- p + theme(axis.text = element_text(size = 18, color = "black"),
               axis.title = element_text(size = 14),
               legend.position = "bottom", legend.box = "horizontal")
p <- p + ylim(c(-4, 4)) + xlim(c(-4, 4))
# file2write <- paste0(dir_out, "scatterplot_snATAC_snRNA_FC." , "1",".png")
# png(file2write, width = 800, height = 900, res = 150)
# print(p)
# dev.off()
file2write <- paste0(dir_out, "scatterplot_snATAC_snRNA_FC.1", ".pdf")
pdf(file2write, width = 5, height = 5, useDingbats = F)
print(p)
dev.off()

# highlight genes with fold change > 1 ----------------------------------------------------
plotdata_df <- plotdata_df %>%
  mutate(highlight = ((avg_log2FC.snATAC >= 1 & avg_log2FC.snRNA >= 1) | (avg_log2FC.snATAC <= -1 & avg_log2FC.snRNA <= -1) | (Gene %in% c("VIM", "FN1", "CDH2", "WNT5B"))))
p <- ggscatter(data = plotdata_df, x = "avg_log2FC.snATAC", y = "avg_log2FC.snRNA", alpha = 0.8, shape = 16, size = 2.5,
               add = "reg.line", add.params = list(color = "grey50", linetype = 2), # Customize reg. line
               conf.int = F # Add confidence interval
)
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey")
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
p <- p + stat_cor(method = "pearson", label.x = 1, label.y = 1.5, size = 4)
p <- p + geom_text_repel(data = subset(x = plotdata_df, highlight == T & avg_log2FC.snATAC > 0),
                         mapping = aes(x = avg_log2FC.snATAC, y = avg_log2FC.snRNA, label = Gene),
                         max.overlaps = Inf, size = 5,
                         segment.size = 0.4, segment.alpha = 0.6, min.segment.length = 0, box.padding = 0.5, color = "red")
p <- p + geom_text_repel(data = subset(x = plotdata_df, highlight == T & avg_log2FC.snATAC < 0),
                         mapping = aes(x = avg_log2FC.snATAC, y = avg_log2FC.snRNA, label = Gene),
                         max.overlaps = Inf, size = 5,
                         segment.size = 0.4, segment.alpha = 0.6, min.segment.length = 0, box.padding = 0.5, color = "blue", xlim = c(-4, -1))
p <- p + guides(color = guide_legend(nrow = 1, override.aes = aes(size = 3), label.theme = element_text(size = 14)))
p <- p + theme_classic()
p <- p + theme(axis.text = element_text(size = 18, color = "black"),
               axis.title = element_text(size = 14),
               legend.position = "bottom", legend.box = "horizontal")
p <- p + ylim(c(-4, 4)) + xlim(c(-4, 4))
p <- p + xlab("Log2(fold change of gene promoter accessiblity)")
p <- p + ylab("Log2(fold change of gene expression)")
# file2write <- paste0(dir_out, "scatterplot_snATAC_snRNA_FC." , "1",".png")
# png(file2write, width = 800, height = 900, res = 150)
# print(p)
# dev.off()
file2write <- paste0(dir_out, "scatterplot_snATAC_snRNA_FC1.wManualAdded", ".pdf")
pdf(file2write, width = 5, height = 5, useDingbats = F)
print(p)
dev.off()

