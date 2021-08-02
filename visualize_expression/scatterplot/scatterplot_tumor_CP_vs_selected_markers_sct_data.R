# Yige Wu @WashU Apr 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA//"
setwd(dir_base)
source("./ccRCC_snRNA_analysis//load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/plotting.R")
source("./ccRCC_snRNA_analysis/plotting_variables.R")
library(ggpubr)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# set parameters ----------------------------------------------------------
gene_x <- "CP"
gene_y <- "VEGFA"
celltype_x <- "Tumor.cells"
celltype_y <- "Tumor.cells"
genes_plot <- c(gene_x, gene_y)

# input dependencies ------------------------------------------------------
## input the average expression
exp_wide_df <- fread(data.table = F, input = "./Resources/Analysis_Results/average_expression/avgexp_sct_data_bycelltype13_bysample_katmai/20210701.v1/33_aliquot_merged.avgexp.SCT.data.bycelltype13_bysample.20210701.v1.tsv")
## input meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210423.v1/meta_data.20210423.v1.tsv")


# make plot data ----------------------------------------------------------
plotdata_wide_df <- exp_wide_df %>%
  filter(V1 %in% genes_plot)
plotdata_raw_df <- melt(data = plotdata_wide_df)
plotdata_raw_df <- plotdata_raw_df %>%
  mutate(id_sample_cell_group = gsub(x = variable, pattern = "SCT\\.", replacement = "")) %>%
  mutate(aliquot = str_split_fixed(string = id_sample_cell_group, pattern = "_", n = 2)[,1]) %>%
  mutate(cell_group.columnname = str_split_fixed(string = id_sample_cell_group, pattern = "_", n = 2)[,2])
plotdata_x_df <- plotdata_raw_df %>%
  filter(V1 == gene_x) %>%
  filter(cell_group.columnname == celltype_x)
plotdata_y_df <- plotdata_raw_df %>%
  filter(V1 == gene_y) %>%
  filter(cell_group.columnname == celltype_y)
plotdata_df <- merge(x = plotdata_x_df %>%
                       select(aliquot, value),
                     y = plotdata_y_df %>%
                       select(aliquot, value),
                     suffixes = c(".x", ".y"), by = c("aliquot"))


# plot --------------------------------------------------------------------
p <- ggplot()
p <- ggscatter(data = plotdata_df, x = "value.x", y = "value.y",
               add = "reg.line",  # Add regressin line
               add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
               conf.int = TRUE # Add confidence interval
)
# Add correlation coefficient
p <- p + stat_cor(method = "spearman", 
                  label.x = quantile(x = plotdata_df$value.x, probs = 0.75, na.rm = T), 
                  label.y = quantile(x = plotdata_df$value.y, probs = 0.75, na.rm = T))
p <- p + xlab(paste0(gene_x, " expression in ", celltype_x))
p <- p + ylab(paste0(gene_y, " expression in ", celltype_y))
p <- p + theme(legend.position = "bottom")
p <- p + guides(color = guide_legend(title = NULL, ncol = 4, label.theme = element_text(size = 8)))
file2write <- paste0(dir_out, celltype_x, "_", gene_x, ".", 
                     celltype_y, "_", gene_y, ".png")
png(file2write, width = 800, height = 800, res = 150)
print(p)
dev.off()
