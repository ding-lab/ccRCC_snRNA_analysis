# Yige Wu @ WashU 2022 May

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "ggplot2"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    install.packages(pkg_name_tmp, dependencies = T)
  }
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg_name_tmp)
  }
  library(package = pkg_name_tmp, character.only = T)
}

# input -------------------------------------------------------------------
enricher_out_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/rna/knockout_cell_lines/pathway_analysis/prioritize_ora_msigdb_h_cp_results_2caki1_cp_vs_2caki1_nt/20220524.v1/ora_msigdb_h_cp_2caki1_cp_vs_2caki1_nt_degs.processed.sig.20220524.v1.tsv")

# make plot data ----------------------------------------------------------
library(dplyr)
plotdata_df <- enricher_out_df %>%
  filter(Keep) %>%
  filter(test == "Caki1_cp_vs_nt.up") %>%
  mutate(size_plot = Count) %>%
  mutate(x_plot = GeneRatio_num*100) %>%
  mutate(log10FDR = -log10(p.adjust)) %>%
  mutate(y_plot = ID)
# pathway_label_df <- data.frame(Description = pathways_selected,
#                                pathway_label = c("SLC transmembrane trasnport", "Eicosanoid metabolism", "Extracellular matrix proteins", "Zinc homeostasis", "Late response to estrogen"))
# 
# plotdata_df$y_plot <- mapvalues(x = plotdata_df$Description, from = pathway_label_df$Description, to = as.vector(pathway_label_df$pathway_label))
plotdata_df <- plotdata_df %>%
  arrange(x_plot)
plotdata_df$y_plot <- factor(x = plotdata_df$y_plot, levels = plotdata_df$y_plot)

# plot enrichment map -----------------------------------------------------
p <- ggplot()
p <- p + geom_point(data = plotdata_df, mapping = aes(x = x_plot, y = y_plot, size = size_plot, color = log10FDR))
# p <- p + scale_color_gradientn(colors = c("blue", "purple", "red"))
p <- p + scale_color_gradientn(colors = c("blue", "purple", "red"),
                               breaks = c(6, 8, 10),
                               guide = guide_colourbar(title = "-log10(p.adjust)", direction = "horizontal", title.position = "top"))
p <- p + scale_size_continuous(breaks = c(40, 60, 80),
                               guide = guide_legend(direction = "horizontal", title = "Gene count", nrow = 1, byrow = T, title.position = "top"))
p <- p + theme_light(base_size = 12)
p <- p + xlab(label = "Gene ratio (%)")
p <- p + xlim(c(0.015, 0.07)*100)
# p <- p + guides(size = guide_legend(title = "Gene count",
#                                     row = 2,
#                                     title.position = "top"))
p <- p + theme(axis.text = element_text(color = "black"),
               axis.title.y = element_blank(), axis.title.x = element_text(size = 10),
               legend.position = "right", legend.box = "vertical")
p


# write output ------------------------------------------------------------
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis//functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
file2write <- paste(dir_out, "dotplot.pdf")
pdf(file2write, width = 7.25, height = 1.5, useDingbats = F)
print(p)
dev.off()
