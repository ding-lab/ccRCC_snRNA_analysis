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
results_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/30ccRCC_tumorcellreclustered/findmarker_30ccRCC_tumorcellreclustered_byres_pairwisebycluster/res.2.tumorcellsreclustered.pairwisebycluster.markers.logfcthreshold.0.25.minpct.0.1.mindiffpct.0.tsv")

# make matrix data -----------------------------------------------------------------
plotdata_df <- results_df %>%
  filter(p_val_adj < 0.05) %>%
  mutate(diff_pct = (pct.1 - pct.2)) %>%
  # filter(diff_pct >= 0.1) %>%
  filter(diff_pct >= 0) %>%
  mutate(x_plot = ifelse(ident.1 > ident.2, ident.1, ident.2)) %>%
  mutate(y_plot = ifelse(ident.1 > ident.2, ident.2, ident.1)) %>%
  group_by(x_plot, y_plot) %>%
  summarise(number_degs = n())

# plotdata_df <- rbind(plotdata_df,
#                      data.frame(x_plot = c(1, 0, 1, 6), y_plot = c(1, 0, 0, 0), number_degs = c(NA, NA, 0, 0)))
plotdata_df$x_plot <- factor(plotdata_df$x_plot)
plotdata_df$y_plot <- factor(plotdata_df$y_plot)
# plotdata_df$number_degs_brks <- cut(plotdata_df$number_degs, breaks = c(-1, 40, 100, 200, 300), labels = c("0-40", "41-100", "101-200", "201-300"))
# plotdata_df$number_degs_brks <- cut(plotdata_df$number_degs, breaks = c(-1, 50, 100*(1:ceiling(max(plotdata_df$number_degs/100)))), labels = c("0-50", "51-100", "101-200", "201-300"))

# plot --------------------------------------------------------------------
p <- ggplot(data = plotdata_df, mapping = aes(x = x_plot, y = y_plot))
p <- p + geom_tile(mapping = aes(fill = number_degs))
p <- p + geom_text(mapping = aes(label = number_degs))
p <- p + scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                               midpoint = 50, space = "Lab",
                               name="Number of\nDEGs")
# p <- p + scale_fill_viridis_d(name="Number of\nDEGs")
p <- p + xlab("Cluster #1") + ylab("Cluster #2")
p <- p + theme_minimal()

# save outputs ------------------------------------------------------------
file2write <- paste0(dir_out, "Number_of_DEGs.png")
png(file2write, width = 1500, height = 1300, res = 150)
print(p)
dev.off()

# 
# # plot --------------------------------------------------------------------
# p <- ggplot(data = plotdata_df, mapping = aes(x = x_plot, y = y_plot))
# p <- p + geom_tile(mapping = aes(fill = number_degs_brks))
# p <- p + geom_text(mapping = aes(label = number_degs), color = "orange")
# # p <- p + scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
# #                                midpoint = 50, limit = c(0, 400), space = "Lab", 
# #                                name="Number of\nDEGs")
# p <- p + scale_fill_viridis_d(name="Number of\nDEGs")
# p <- p + xlab("Cluster #1") + ylab("Cluster #2")
# p <- p + theme_minimal()
# 
# # save outputs ------------------------------------------------------------
# file2write <- paste0(dir_out, "Number_of_DEGs.png")
# png(file2write, width = 1000, height = 800, res = 150)
# print(p)
# dev.off()
