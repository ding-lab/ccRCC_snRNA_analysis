# Yige Wu @WashU Oct 2020

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
## input clinical data
specimen_clinical_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/extract_specimen_clinical_data/20200717.v1/snRNA_ccRCC_Specimen_Clinicl_Data.20200717.v1.tsv")
## input cell type fraction
cellgroupfrac_df <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/count_fraction/count_celltypeshorter_fraction_per_sample/20201203.v1/CellGroupBarcodes_Number_and_Fraction_per_Sample20201203.v1.tsv")


# set plotting parameters -------------------------------------------------
symnum.args <- list(cutpoints = c(0, 0.0001, 0.001, 0.01, 0.05, 1), symbols = c("****", "***", "**", "*", "ns"))
pos <- position_jitter(width = 0.2, seed = 1)
my_comparisons <- list(c("G1/2", "G3"),c("G3", "G4"),c("G1/2", "G4"))

# plot by cell group ------------------------------------------------------
table(cellgroupfrac_df$Cell_group)
cellgroup_tmp <- "CD8 CTL exhausted"
cellgroup_tmp <- "Mast cells"

for (cellgroup_tmp in unique(cellgroupfrac_df$Cell_group)) {
  ## make plot data
  plot_data_df <- specimen_clinical_df %>%
    filter(Sample_Type == "Tumor") %>%
    filter(snRNA_available == TRUE) %>%
    select(Aliquot.snRNA.WU, Histologic_Grade)
  plot_data_df <- merge(x = plot_data_df,
                        y = cellgroupfrac_df %>%
                          filter(Cell_group == cellgroup_tmp) %>%
                          select(Aliquot_WU, Frac_CellGroupBarcodes_ByAliquot),
                        by.x = c("Aliquot.snRNA.WU"), by.y = c("Aliquot_WU"), all.x = T)
  plot_data_df <- plot_data_df %>%
    mutate(x_plot = ifelse(Histologic_Grade %in% c("G1", "G2"), "G1/2", Histologic_Grade)) %>%
    mutate(y_plot = ifelse(is.na(Frac_CellGroupBarcodes_ByAliquot), 0, Frac_CellGroupBarcodes_ByAliquot))
  ## plot
  p = ggplot(plot_data_df, aes(x=x_plot, y=y_plot))
  p = p + geom_violin(aes(fill = x_plot),  color = NA, alpha = 1)
  p = p + geom_boxplot(width=.1)
  p = p + geom_point(color = "black", fill = "black", shape = 16, position = pos, stroke = 0, alpha = 0.6, size = 1)
  p = p + stat_compare_means(data = plot_data_df, 
                             mapping = aes(x = x_plot, y = y_plot, label = ..p.signif..), comparisons = my_comparisons,
                             symnum.args = symnum.args,
                             hide.ns = F, method = "wilcox.test")
  p <- p + theme(legend.position = "none")
  p <- p + theme(axis.title = element_blank())
  file2write <- paste0(dir_out, gsub(x = cellgroup_tmp, pattern = "\\/", replacement = "_"), ".png")
  png(file2write, width = 600, height = 600, res = 150)
  print(p)
  dev.off()
}




