# Yige Wu @WashU Apr 2022
## plot cell type on integration UMAP

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "ggplot2"
)
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input DEGs
results_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/30ccRCC_tumorcellreclustered/findallmarker_30ccRCC_tumorcellreclustered_byres/20220408.v1/tumorcellsreclustered.markers.logfcthreshold.0.25.minpct.0.1.mindiffpct.0.byresolution.tsv")

# make matrix data -----------------------------------------------------------------
plotdata_df <- results_df %>%
  filter(p_val_adj < 0.05) %>%
  mutate(diff_pct = (pct.1 - pct.2)) %>%
  # filter(diff_pct >= 0.1) %>%
  filter(diff_pct >= 0) %>%
  group_by(resolution, cluster) %>%
  summarise(number_degs = n())
plotdata_df$x_plot <- factor(plotdata_df$resolution)
textdata_df <- plotdata_df %>%
  group_by(resolution) %>%
  summarise(number_clusters = n()) %>%
  mutate(text = paste0("n = ", number_clusters))
textdata_df$x_plot <- factor(textdata_df$resolution)
# plot --------------------------------------------------------------------
pos <- position_jitter(width = 0.1, seed = 1)
p <- ggplot(data = plotdata_df, mapping = aes(x = x_plot, y = number_degs))
p <- p + stat_boxplot(geom = "errorbar", width = 0.2)
p <- p + geom_boxplot(mapping = aes(group = x_plot), outlier.color = "red")
p <- p + geom_point(position = pos, alpha = 0.8)
p <- p + geom_hline(yintercept = 20, linetype = 2)
p <- p + geom_text(data = textdata_df, mapping = aes(x = x_plot, y = 410, label = text), size = 6.5)
p <- p + scale_y_continuous(breaks = c(0, 20, 100, 200, 300, 400))
p <- p + xlab("Resolution") + ylab("Number of unique markers per cluster")
p <- p + theme_classic(base_size = 15)
p <- p + theme(axis.text = element_text(color = "black", size = 20), axis.title = element_text(color = "black", size = 20))

# save outputs ------------------------------------------------------------
file2write <- paste0(dir_out, "Number_of_DEGs_byresolution.pdf")
pdf(file2write, width = 6, height = 6, useDingbats = F)
print(p)
dev.off()

file2write <- paste0(dir_out, "Number_of_DEGs_byresolution.png")
png(file2write, width = 1000, height = 800, res = 150)
print(p)
dev.off()
