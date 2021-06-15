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

## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
peaks2degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/bap1/overlap_bap1_associated_peaks_with_degs/20210611.v1/BAP1_DEG_associated_Peaks.20210611.v1.tsv")

# make plot data ----------------------------------------------------------
plotdata_df <- peaks2degs_df %>%
  filter(!is.na(avg_log2FC.snATAC) & !is.na(avg_log2FC.snRNA))

# plot --------------------------------------------------------------------
p <- ggscatter(data = plotdata_df, x = "avg_log2FC.snATAC", y = "avg_log2FC.snRNA", color = "DAP_type",
                   add = "reg.line",  # Add regressin line
                   add.params = list(color = "blue", fill = "lightgray", linetype = 2), # Customize reg. line
                   conf.int = F # Add confidence interval
)
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey")
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
p <- p + stat_cor(method = "pearson", label.x = 0.5, label.y = 1.2)
p <- p + guides(color = guide_legend(nrow = 2, override.aes = aes(size = 3), label.theme = element_text(size = 14)))
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

