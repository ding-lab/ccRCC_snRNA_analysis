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
library(ggrastr)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input merged tumor content from bulk and snRNA
celltype_frac_snrna_df <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/count_fraction/count_celltype_wepithelial_fraction_per_sample/20211007.v1/CellGroupBarcodes_Number_and_Fraction_per_Sample20211007.v1.tsv")
celltype_frac_snatac_df <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/count_fraction/count_celltype_wepithelial_fraction_per_sample_snatac/20211007.v1/CellGroupBarcodes_Number_and_Fraction_per_snATAc_Sample20211007.v1.tsv")

# make plot data ----------------------------------------------------------
plot_data_df <- merge(x = celltype_frac_snrna_df, y = celltype_frac_snatac_df, by = c("Aliquot_WU", "Cell_group"), all.x = T, suffixes = c(".snRNA", ".snATAC"))
plot_data_df <- plot_data_df %>%
  mutate(x_plot = Frac_CellGroupBarcodes_ByAliquot.snRNA) %>%
  mutate(y_plot = Frac_CellGroupBarcodes_ByAliquot.snATAC) %>%
  filter(!is.na(x_plot) & !is.na(y_plot))

# make scatterplot --------------------------------------------------------
## reference: https://rpkgs.datanovia.com/ggpubr/reference/stat_cor.html
p <- ggplot()
p <- p + geom_point_rast(data = plot_data_df, mapping = aes(x = x_plot, y = y_plot), alpha = 0.8, shape = 16)
# p <- ggscatter(plot_data_df, x = "x_plot", y = "y_plot", alpha = 0.8,
#                add = "reg.line",  # Add regressin line
#                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
#                conf.int = TRUE # Add confidence interval
# )
# Add correlation coefficient
p <- p + stat_cor(method = "pearson", label.x = 0.01, label.y = 0.4, size = 6)
p <- p + xlab("Cell type content estimated from snRNA Data")
p <- p + ylab("Cell type content estimated from snATAC Data")
p <- p + theme_classic(base_size = 16)
p <- p + theme(axis.text = element_text(color = "black", size = 16),
               axis.title = element_text(color = "black", size = 16))
p

# save scatterplot --------------------------------------------------------
file2write <- paste0(dir_out, "scatterplot",".png")
png(file2write, width = 800, height = 800, res = 150)
print(p)
dev.off()
file2write <- paste0(dir_out, "scatterplot",  ".pdf")
pdf(file2write, width = 5, height = 5, useDingbats = F)
print(p)
dev.off()


