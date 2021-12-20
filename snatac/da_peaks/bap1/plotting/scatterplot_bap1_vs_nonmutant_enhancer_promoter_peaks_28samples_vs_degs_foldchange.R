# Yige Wu @WashU Jun 2021
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
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
peaks2degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/bap1/overlap_degs/overlap_bap1_vs_nonmutant_enhancer_promoter_peaks_28samples_with_degs/20211011.v1/BAP1_vs_NonMutant_DAP_DEG_Merged.20211011.v1.tsv")

# make plot data ----------------------------------------------------------
plotdata_df <- peaks2degs_df %>%
  filter(!is.na(avg_log2FC.snATAC) & !is.na(avg_log2FC.snRNA)) %>%
  select(avg_log2FC.snATAC, avg_log2FC.snRNA, Gene, peak2gene_type) %>%
  unique()
## make colors
colors_peak2genetype <- brewer.pal(n = 7, name = "Dark2")[c(4, 6)]
names(colors_peak2genetype) <- c("Promoter", "Enhancer")

# plot highlight genes ----------------------------------------------------
plotdata_df <- plotdata_df %>%
  mutate(highlight = (Gene %in% c("PTPRJ", "DLC1", "DDIT4", "PEBP1")  | (Gene %in% c("SLC38A1", "RAPGEF5", "EPHA6", "DUSP1", "FABP6") & avg_log2FC.snATAC > 0) | (Gene == "CES3" & avg_log2FC.snATAC < -1.5)))

p <- ggscatter(data = plotdata_df, x = "avg_log2FC.snATAC", y = "avg_log2FC.snRNA", color = "peak2gene_type", alpha = 0.8, shape = 16, size = 2.5,
               add = "reg.line", add.params = list(color = "grey50", linetype = 2), # Customize reg. line
               conf.int = F # Add confidence interval
)
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey")
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
# p <- p + stat_cor(method = "pearson", label.x = 0, label.y = 0, size = 4)
p <- p + geom_text_repel(data = subset(x = plotdata_df, highlight == T), 
                         mapping = aes(x = avg_log2FC.snATAC, y = avg_log2FC.snRNA, label = Gene), 
                         max.overlaps = Inf, size = 6,
                         segment.size = 0.4, segment.alpha = 0.6, min.segment.length = 0, box.padding = 0.5)
p <-  p + scale_color_manual(values = colors_peak2genetype)
p <- p + guides(color = guide_legend(nrow = 1, override.aes = aes(size = 3), label.theme = element_text(size = 14)))
p <- p + theme(axis.text = element_text(size = 18, color = "black"),
               axis.title = element_text(size = 14),
               legend.position = "bottom", legend.box = "horizontal")
# file2write <- paste0(dir_out, "scatterplot_snATAC_snRNA_FC." , "1",".png")
# png(file2write, width = 800, height = 900, res = 150)
# print(p)
# dev.off()
file2write <- paste0(dir_out, "scatterplot_snATAC_snRNA_FC.","PTPRJ_CES3",".pdf")
pdf(file2write, width = 5.5, height = 6, useDingbats = F)
print(p)
dev.off()

# # # plot --------------------------------------------------------------------
# p <- ggscatter(data = plotdata_df, x = "avg_log2FC.snATAC", y = "avg_log2FC.snRNA", color = "peak2gene_type",
#                    add = "reg.line",  # Add regressin line
#                    add.params = list(color = "blue", fill = "lightgray", linetype = 2), # Customize reg. line
#                    conf.int = F # Add confidence interval
# )
# p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey")
# p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
# p <- p + stat_cor(method = "pearson", label.x = 0, label.y = 0)
# p <- p + guides(color = guide_legend(nrow = 1, override.aes = aes(size = 3), label.theme = element_text(size = 14)))
# p <- p + theme(axis.text = element_text(size = 14),
#                axis.title = element_text(size = 14),
#                legend.position = "bottom", legend.box = "horizontal")
# # file2write <- paste0(dir_out, "scatterplot_snATAC_snRNA_FC.",".png")
# # png(file2write, width = 800, height = 900, res = 150)
# # print(p)
# # dev.off()
# file2write <- paste0(dir_out, "scatterplot_snATAC_snRNA_FC.",".pdf")
# pdf(file2write, width = 5.5, height = 6, useDingbats = F)
# print(p)
# dev.off()
# 
# # plot highlight genes ----------------------------------------------------
# plotdata_df <- plotdata_df %>%
#   mutate(highlight = (avg_log2FC.snATAC*avg_log2FC.snRNA >0) & (abs(avg_log2FC.snATAC) > 1 | abs(avg_log2FC.snRNA) > 1))
# 
# p <- ggscatter(data = plotdata_df, x = "avg_log2FC.snATAC", y = "avg_log2FC.snRNA", color = "peak2gene_type",
#                add = "reg.line",  # Add regressin line
#                add.params = list(color = "grey50", linetype = 2), # Customize reg. line
#                conf.int = F # Add confidence interval
# )
# p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey")
# p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
# # p <- p + stat_cor(method = "pearson", label.x = 0, label.y = 0)
# p <- p + geom_text_repel(data = subset(x = plotdata_df, highlight == T), 
#                          mapping = aes(x = avg_log2FC.snATAC, y = avg_log2FC.snRNA, label = Gene), 
#                          max.overlaps = Inf, min.segment.length = 0)
# p <- p + guides(color = guide_legend(nrow = 2))
# p <- p + theme(axis.text = element_text(size = 14),
#                axis.title = element_text(size = 14),
#                legend.position = "bottom", legend.box = "horizontal")
# file2write <- paste0(dir_out, "scatterplot_snATAC_snRNA_FC." , "1",".png")
# # png(file2write, width = 800, height = 900, res = 150)
# # print(p)
# # dev.off()
# file2write <- paste0(dir_out, "scatterplot_snATAC_snRNA_FC.","1",".pdf")
# pdf(file2write, width = 5.5, height = 6, useDingbats = F)
# print(p)
# dev.off()
# 
# # plot highlight genes ----------------------------------------------------
# plotdata_df <- plotdata_df %>%
#   mutate(highlight = (avg_log2FC.snATAC*avg_log2FC.snRNA > 0) & (abs(avg_log2FC.snATAC) > 0.5 & abs(avg_log2FC.snRNA) > 0.5))
# 
# p <- ggscatter(data = plotdata_df, x = "avg_log2FC.snATAC", y = "avg_log2FC.snRNA", color = "peak2gene_type",
#                add = "reg.line",  # Add regressin line
#                add.params = list(color = "grey50", linetype = 2), # Customize reg. line
#                conf.int = F # Add confidence interval
# )
# p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey")
# p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
# # p <- p + stat_cor(method = "pearson", label.x = 0, label.y = 0)
# p <- p + geom_text_repel(data = subset(x = plotdata_df, highlight == T), 
#                          mapping = aes(x = avg_log2FC.snATAC, y = avg_log2FC.snRNA, label = Gene), 
#                          max.overlaps = Inf, min.segment.length = 0)
# p <- p + guides(color = guide_legend(nrow = 2))
# p <- p + theme(axis.text = element_text(size = 14),
#                axis.title = element_text(size = 14),
#                legend.position = "bottom", legend.box = "horizontal")
# # file2write <- paste0(dir_out, "scatterplot_snATAC_snRNA_FC." , "1",".png")
# # png(file2write, width = 800, height = 900, res = 150)
# # print(p)
# # dev.off()
# file2write <- paste0(dir_out, "scatterplot_snATAC_snRNA_FC.","0.5",".pdf")
# pdf(file2write, width = 5.5, height = 6, useDingbats = F)
# print(p)
# dev.off()

# # plot highlight genes ----------------------------------------------------
# plotdata_df <- plotdata_df %>%
#   mutate(highlight = (Gene %in% c("PTPRJ", "CES3")))
# 
# p <- ggscatter(data = plotdata_df, x = "avg_log2FC.snATAC", y = "avg_log2FC.snRNA", color = "peak2gene_type",
#                add = "reg.line",  # Add regressin line
#                add.params = list(color = "grey50", linetype = 2, alpha = 0.8), # Customize reg. line
#                conf.int = F # Add confidence interval
# )
# p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey")
# p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
# p <- p + stat_cor(method = "pearson", label.x = 0, label.y = 0, size = 4)
# p <- p + geom_text_repel(data = subset(x = plotdata_df, highlight == T), 
#                          mapping = aes(x = avg_log2FC.snATAC, y = avg_log2FC.snRNA, label = Gene), 
#                          max.overlaps = Inf, size = 6,
#                          segment.size = 0.4, segment.alpha = 0.6, min.segment.length = 0)
# p <- p + guides(color = guide_legend(nrow = 1, override.aes = aes(size = 3), label.theme = element_text(size = 14)))
# p <- p + theme(axis.text = element_text(size = 14),
#                axis.title = element_text(size = 14),
#                legend.position = "bottom", legend.box = "horizontal")
# # file2write <- paste0(dir_out, "scatterplot_snATAC_snRNA_FC." , "1",".png")
# # png(file2write, width = 800, height = 900, res = 150)
# # print(p)
# # dev.off()
# file2write <- paste0(dir_out, "scatterplot_snATAC_snRNA_FC.","PTPRJ_CES3",".pdf")
# pdf(file2write, width = 5.5, height = 6, useDingbats = F)
# print(p)
# dev.off()

