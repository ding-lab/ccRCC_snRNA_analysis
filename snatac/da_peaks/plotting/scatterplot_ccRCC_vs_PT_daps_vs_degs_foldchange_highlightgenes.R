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
peaks2degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/overlap_deg/overlap_ccRCC_vs_PT_DA_peaks_with_degs/20210614.v1/ccRCC_vs_PT_DEG_associated_DA_Peaks.20210614.v1.tsv")

# make plot data ----------------------------------------------------------
plotdata_df <- peaks2degs_df %>%
  filter(!is.na(avg_log2FC.snATAC) & !is.na(avg_log2FC.snRNA))

# plot --------------------------------------------------------------------
p <- ggplot()
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey")
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
p <- p + geom_vline(xintercept = 1, linetype = 2, color = "red")
p <- p + geom_hline(yintercept = 1, linetype = 2, color = "red")
p <- p + geom_point(data = plotdata_df, mapping = aes(x = avg_log2FC.snATAC, y = avg_log2FC.snRNA, color = DAP_type))
p <- p + geom_text_repel(data = subset(x = plotdata_df, abs(avg_log2FC.snATAC) > 1 & abs(avg_log2FC.snRNA) > 1), 
                         mapping = aes(x = avg_log2FC.snATAC, y = avg_log2FC.snRNA, label = Gene), 
                         max.overlaps = Inf, min.segment.length = 0)
p <- p + theme_classic()
p <- p + theme(legend.position = "top")
p <- p + guides(color = guide_legend(nrow = 2))
file2write <- paste0(dir_out, "scatterplot_snATAC_snRNA_FC.",".png")
png(file2write, width = 800, height = 900, res = 150)
print(p)
dev.off()

