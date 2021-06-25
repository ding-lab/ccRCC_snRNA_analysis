# Yige Wu @WashU Jun 2021

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
peaks2degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/pbrm1/overlap_degs/overlap_pbrm1_vs_nonmutant_enhancer_promoter_peaks_with_degs/20210625.v1/PBRM1_vs_NonMutant_DAP2DEG.20210625.v1.tsv")

# make plot data ----------------------------------------------------------
plotdata_df <- peaks2degs_df %>%
  filter(!is.na(avg_log2FC.snATAC) & !is.na(avg_log2FC.snRNA)) %>%
  select(avg_log2FC.snATAC, avg_log2FC.snRNA, Gene, peak2gene_type) %>%
  unique()

# plot --------------------------------------------------------------------
p <- ggscatter(data = plotdata_df, x = "avg_log2FC.snATAC", y = "avg_log2FC.snRNA", color = "peak2gene_type",
                   add = "reg.line",  # Add regressin line
                   add.params = list(color = "grey50",linetype = 2), # Customize reg. line
                   conf.int = F # Add confidence interval
)
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey")
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
p <- p + stat_cor(method = "pearson", label.x = 0, label.y = 0)
p <- p + guides(color = guide_legend(nrow = 1, override.aes = aes(size = 3), label.theme = element_text(size = 14)))
p <- p + theme(axis.text = element_text(size = 14),
               axis.title = element_text(size = 14),
               legend.position = "bottom", legend.box = "horizontal")
# file2write <- paste0(dir_out, "scatterplot_snATAC_snRNA_FC.",".png")
# png(file2write, width = 800, height = 900, res = 150)
# print(p)
# dev.off()
file2write <- paste0(dir_out, "scatterplot_snATAC_snRNA_FC.",".pdf")
pdf(file2write, width = 5.5, height = 6, useDingbats = F)
print(p)
dev.off()

# plot highlight genes ----------------------------------------------------
plotdata_df <- plotdata_df %>%
  mutate(highlight = (avg_log2FC.snATAC*avg_log2FC.snRNA >0) & (abs(avg_log2FC.snATAC) > 1 | abs(avg_log2FC.snRNA) > 1))

p <- ggscatter(data = plotdata_df, x = "avg_log2FC.snATAC", y = "avg_log2FC.snRNA", color = "peak2gene_type",
               add = "reg.line",  # Add regressin line
               add.params = list(color = "grey50", linetype = 2), # Customize reg. line
               conf.int = F # Add confidence interval
)
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey")
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
p <- p + stat_cor(method = "pearson", label.x = 0, label.y = 0)
p <- p + geom_text_repel(data = subset(x = plotdata_df, highlight == T), 
                         mapping = aes(x = avg_log2FC.snATAC, y = avg_log2FC.snRNA, label = Gene), 
                         max.overlaps = Inf, min.segment.length = 0)
p <- p + guides(color = guide_legend(nrow = 2))
p <- p + theme(axis.text = element_text(size = 14),
               axis.title = element_text(size = 14),
               legend.position = "bottom", legend.box = "horizontal")
# file2write <- paste0(dir_out, "scatterplot_snATAC_snRNA_FC." , "1",".png")
# png(file2write, width = 800, height = 900, res = 150)
# print(p)
# dev.off()
file2write <- paste0(dir_out, "scatterplot_snATAC_snRNA_FC.","1",".pdf")
pdf(file2write, width = 5.5, height = 6, useDingbats = F)
print(p)
dev.off()

# plot highlight genes ----------------------------------------------------
plotdata_df <- plotdata_df %>%
  mutate(highlight = T)

p <- ggscatter(data = plotdata_df, x = "avg_log2FC.snATAC", y = "avg_log2FC.snRNA", color = "peak2gene_type",
               add = "reg.line",  # Add regressin line
               add.params = list(color = "grey50", linetype = 2), # Customize reg. line
               conf.int = F # Add confidence interval
)
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey")
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
# p <- p + stat_cor(method = "pearson", label.x = 0, label.y = 0)
p <- p + geom_text_repel(data = subset(x = plotdata_df, highlight == T), 
                         mapping = aes(x = avg_log2FC.snATAC, y = avg_log2FC.snRNA, label = Gene), 
                         max.overlaps = Inf, min.segment.length = 0)
p <- p + guides(color = guide_legend(nrow = 2))
p <- p + theme(axis.text = element_text(size = 14),
               axis.title = element_text(size = 14),
               legend.position = "bottom", legend.box = "horizontal")
# file2write <- paste0(dir_out, "scatterplot_snATAC_snRNA_FC." , "1",".png")
# png(file2write, width = 800, height = 900, res = 150)
# print(p)
# dev.off()
file2write <- paste0(dir_out, "scatterplot_snATAC_snRNA_FC.","highlightall",".pdf")
pdf(file2write, width = 5.5, height = 6, useDingbats = F)
print(p)
dev.off()

