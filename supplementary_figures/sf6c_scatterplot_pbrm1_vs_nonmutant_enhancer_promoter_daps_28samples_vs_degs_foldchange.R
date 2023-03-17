# Yige Wu @WashU Jun 2021

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
  "circlize"
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

# input dependencies ------------------------------------------------------
peaks2degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/pbrm1/overlap_degs/overlap_pbrm1_vs_nonmutant_enhancer_promoter_peaks_28samples_with_degs/20211011.v1/PBRM1_vs_NonMutant_DAP2DEG.20211011.v1.tsv")

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
  mutate(highlight = (Gene %in% c("PBX1", "RP1", "NCALD")) | (Gene %in% c("ANGPTL4", "CAV1", "PPFIA4", "CAV1", "CAV2", "TNFAIP8", "EDN1")))

p <- ggscatter(data = plotdata_df, x = "avg_log2FC.snATAC", y = "avg_log2FC.snRNA", color = "peak2gene_type",
               add = "reg.line",  # Add regressin line
               add.params = list(color = "grey50", linetype = 2), # Customize reg. line
               conf.int = F # Add confidence interval
)
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey")
p <- p + geom_hline(yintercept = 0, linetype = 2, color = "grey")
p <- p + stat_cor(method = "pearson", label.x = -1, label.y = 0)
p <-  p + scale_color_manual(values = colors_peak2genetype)
p <- p + geom_text_repel(data = subset(x = plotdata_df, highlight == T), 
                         mapping = aes(x = avg_log2FC.snATAC, y = avg_log2FC.snRNA, label = Gene), 
                         max.overlaps = Inf, size = 6, 
                         segment.size = 0.4, segment.alpha = 0.6, min.segment.length = 0, box.padding = 0.5)
p <- p + xlab("Log2(fold change of peak accessibility)")
p <- p + ylab("Log2(fold change of gene expression)")
p <- p + theme_classic()
p <- p + guides(color = guide_legend(nrow = 1))
p <- p + theme(axis.text = element_text(size = 18, color = "black"),
               axis.title = element_text(size = 18),
               legend.position = "bottom", legend.box = "horizontal")
file2write <- paste0(dir_out, "scatterplot_snATAC_snRNA_FC.","selectedGenes",".pdf")
pdf(file2write, width = 5.5, height = 6, useDingbats = F)
print(p)
dev.off()
## save source data
source_data_df <- plotdata_df %>%
  select(avg_log2FC.snATAC, avg_log2FC.snRNA, peak2gene_type)
write.table(x = source_data_df, file = "~/Desktop/SF6b.SourceData.tsv", quote = F, sep = "\t", row.names = F)

