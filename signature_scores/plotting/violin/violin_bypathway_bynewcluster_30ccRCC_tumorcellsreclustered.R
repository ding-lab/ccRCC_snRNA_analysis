# Yige Wu @WashU Mar 2021
## https://satijalab.org/seurat/archive/v3.0/integration.html
## also used references

# set up libraries and output directory -----------------------------------
## getting the path to the current script
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
path_this_script <- thisFile()
## set working directory
# dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "ggplot2",
  "Polychrome"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)


# input ------------------------------------------------------
## input signature scores by barcode
sigScores <- readRDS("./Resources/Analysis_Results/signature_scores/run_vision/run_vision_on_30ccRCC_tumorcellreclustered/20220406.v1/ccRCC.30ccRCC.TumorCellsReclustered.Vision.sigScores.20220406.v1.RDS")
print("Finish reading the sigScores matrix!\n")
## input the barcode to new cluster
barcode2cluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/integration/seuratintegrate_34_ccRCC_samples/annotate_clusterid_on_30ccRCCtumorcellreclustered_byresolution/20220406.v1/ccRCC.34Sample.Tumorcells.Integrated.ReciprocalPCA.Metadata.20220406.v1.tsv")
## specify pathways to plot
genesets_plot <- c("HALLMARK_FATTY_ACID_METABOLISM", "KEGG_FATTY_ACID_METABOLISM", "WP_FATTY_ACID_BIOSYNTHESIS", 
                   "KEGG_GLYCOLYSIS_GLUCONEOGENESIS", "HALLMARK_GLYCOLYSIS",
                   "HALLMARK_OXIDATIVE_PHOSPHORYLATION", "KEGG_OXIDATIVE_PHOSPHORYLATION", "WP_OXIDATIVE_PHOSPHORYLATION",
                   "GOBP_KIDNEY_MORPHOGENESIS", 
                   "GOBP_ADHERENS_JUNCTION_ASSEMBLY", "KEGG_ADHERENS_JUNCTION",
                   "HALLMARK_ANGIOGENESIS", "WP_ANGIOGENESIS")
genesets_plot <- genesets_plot[genesets_plot %in% colnames(sigScores)]; genesets_plot

# prepare data to plot -----------------------------------------------------------------
scores_wide_df <- data.frame(sigScores[,genesets_plot])
colnames(scores_wide_df) <- genesets_plot
scores_wide_df$barcode <- rownames(sigScores)
scores_long_df <- reshape2::melt(data = scores_wide_df, id.vars = c("barcode"))
scores_long_df$clusterid_new <- mapvalues(x = scores_long_df$barcode, from = barcode2cluster_df$barcode, to = as.vector(barcode2cluster_df$clusterid_new))
head(scores_long_df)
clusterids <- sort(unique(scores_long_df$clusterid_new))
colors_bycluster <- Polychrome::palette36.colors(n = (length(clusterids)+1))[-2]
names(colors_bycluster) <- clusterids

# plot --------------------------------------------------------------
geneset_tmp <- genesets_plot[1]
for (geneset_tmp in genesets_plot) {
  plotdata_df <- scores_long_df %>%
    filter(variable == geneset_tmp) %>%
    mutate(x_plot = clusterid_new) %>%
    mutate(y_plot = value)
  ## plot
  p = ggplot(plotdata_df, aes(x=x_plot, y=y_plot))
  p = p + geom_violin(aes(fill = x_plot),  color = NA, alpha = 1)
  p = p + geom_boxplot(width=.1)
  p <- p + scale_fill_manual(values = colors_bycluster)
  # p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.8, size = 0.1)
  # p = p + stat_compare_means(data = plot_data_df, 
  #                            mapping = aes(x = x_plot, y = y_plot), 
  #                            hide.ns = F, method = "wilcox.test")
  p <- p + ggtitle(label = paste0(geneset_tmp), subtitle = "signature scores by tumor-cell cluster")
  p <- p + theme_classic()
  p <- p + theme(legend.position = "none")
  p <- p + theme(axis.title = element_blank())
  file2write <- paste0(dir_out, geneset_tmp, ".png")
  png(file2write, width = 600, height = 600, res = 150)
  print(p)
  dev.off()
}
