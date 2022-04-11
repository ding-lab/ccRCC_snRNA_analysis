# Yige Wu @WashU May 2020
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
results_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/30ccRCC_tumorcellreclustered/findallmarker_30ccRCC_tumorcellreclustered_byres/res.0.5.tumorcellsreclustered.markers.logfcthreshold.0.25.minpct.0.1.mindiffpct.0.tsv")
results_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/30ccRCC_tumorcellreclustered/findallmarker_30ccRCC_tumorcellreclustered_byres/res.1.tumorcellsreclustered.markers.logfcthreshold.0.25.minpct.0.1.mindiffpct.0.tsv")
results_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/30ccRCC_tumorcellreclustered/findallmarker_30ccRCC_tumorcellreclustered_byres/res.2.tumorcellsreclustered.markers.logfcthreshold.0.25.minpct.0.1.mindiffpct.0.tsv")
results_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/30ccRCC_tumorcellreclustered/findallmarker_30ccRCC_tumorcellreclustered_byres/res.3.tumorcellsreclustered.markers.logfcthreshold.0.25.minpct.0.1.mindiffpct.0.tsv")

# make matrix data -----------------------------------------------------------------
count_df <- as.data.frame(results_df) %>%
  # filter(p_val_adj < 0.05) %>%
  filter(p_val_adj < 0.05) %>%
  mutate(diff_pct = (pct.1 - pct.2)) %>%
  filter(diff_pct >= 0) %>%
  group_by(cluster) %>%
  summarise(number_degs = n())


plotdata_df$x_plot <- factor(plotdata_df$x_plot)
plotdata_df$y_plot <- factor(plotdata_df$y_plot)

# plot --------------------------------------------------------------------
p <- ggplot(data = plotdata_df, mapping = aes(x = x_plot, y = y_plot))
p <- p + geom_tile(mapping = aes(fill = number_degs))
p <- p + geom_text(mapping = aes(label = number_degs))
p <- p + scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                               midpoint = 50, limit = c(0, 400), space = "Lab", 
                               name="Number of\nDEGs")
p <- p + xlab("Cluster #1") + ylab("Cluster #2")
p <- p + theme_minimal()
p
