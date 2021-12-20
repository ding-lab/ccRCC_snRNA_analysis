# Yige Wu @WashU March 2021

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
## input merged tumor content from bulk and snRNA
merged_tumorcontent_df <- fread(input = "./Resources/Analysis_Results/bulk/tumor_content/merge_tumor_content_from_bulk_and_snRNA/20211007.v1/Perc_Tumor_Content_from_snRNA_and_bulkRNA.20211007.v1.tsv", data.table = F)

# make plot data ----------------------------------------------------------
mean(merged_tumorcontent_df$Frac_CellGroupBarcodes_ByAliquot)
plot_data_df <- merged_tumorcontent_df
plot_data_df <- plot_data_df %>%
  mutate(x_plot = Frac_CellGroupBarcodes_ByAliquot) %>%
  mutate(y_plot = ESTIMATE_TumorPurity_RNA) %>%
  filter(!is.na(x_plot) & !is.na(y_plot))

# make scatterplot --------------------------------------------------------
## reference: https://rpkgs.datanovia.com/ggpubr/reference/stat_cor.html
p <- ggplot()
# p <- p + geom_point_rast(data = plot_data_df, mapping = aes(x = x_plot, y = y_plot), alpha = 0.8, shape = 16)
p <- ggscatter(plot_data_df, x = "x_plot", y = "y_plot", alpha = 0.8,
               add = "reg.line",  # Add regressin line
               add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
               conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
p <- p + stat_cor(method = "pearson", label.x = 0.4, label.y = 0.8, size = 6)
p <- p + xlab("Tumor Content Estimated from snRNA Data")
p <- p + ylab("Tumor Content Estimated from Bulk RNA Data")
p <- p + theme_classic(base_size = 16)
p <- p + theme(axis.text = element_text(color = "black", size = 16),
               axis.title = element_text(color = "black", size = 16))
p

# save scatterplot --------------------------------------------------------
file2write <- paste0(dir_out, "scatterplot_tumor_content_from_bulk_and_snRNA.",".png")
png(file2write, width = 800, height = 800, res = 150)
print(p)
dev.off()
file2write <- paste0(dir_out, "scatterplot_tumor_content_from_bulk_and_snRNA",  ".pdf")
pdf(file2write, width = 5, height = 5, useDingbats = F)
print(p)
dev.off()


