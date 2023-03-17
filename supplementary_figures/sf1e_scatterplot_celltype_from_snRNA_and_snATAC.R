# Yige Wu @WashU March 2021

# set up libraries and output directory -----------------------------------
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "ggrastr",
  "ggplot2",
  "ggpubr"
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
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
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
  filter(!is.na(x_plot) & !is.na(y_plot)) %>%
  select(x_plot, y_plot)

# make scatterplot --------------------------------------------------------
## reference: https://rpkgs.datanovia.com/ggpubr/reference/stat_cor.html
p <- ggplot(data = plot_data_df, mapping = aes(x = x_plot, y = y_plot))
p <- p + geom_smooth(method=lm , color="blue", se=FALSE, alpha = 0.5)
p <- p + geom_point_rast(alpha = 0.8, shape = 16)
# Add correlation coefficient
p <- p + stat_cor(method = "pearson", label.x = 0.01, label.y = 0.4, size = 6)
p <- p + xlab("Cell type content estimated from snRNA Data")
p <- p + ylab("Cell type content estimated from snATAC Data")
p <- p + theme_classic(base_size = 16)
p <- p + theme(axis.text = element_text(color = "black", size = 16),
               axis.title = element_text(color = "black", size = 16))

lm_fit <- lm(y_plot ~ x_plot, data = plot_data_df)
summary(lm_fit)
lm_fit$coefficients
summary(lm_fit)$coefficients[,4]  

# save scatterplot --------------------------------------------------------
file2write <- paste0(dir_out, "scatterplot",  ".pdf")
pdf(file2write, width = 5, height = 5, useDingbats = F)
print(p)
dev.off()
write.table(x = plot_data_df, file = "~/Desktop/SF1e.Bottom.SourceData.tsv", quote = F, sep = "\t", row.names = F)


