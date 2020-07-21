# Yige Wu @WashU Jul 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input number of genes per cell ------------------------------------------
metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/individual_sample/fetch_data/fetch_data/20200717.v1/Barcode2MetaData.20200717.v1.tsv")

# make data for plotting --------------------------------------------------
plot_data_df <- metadata_df %>%
  select(nFeature_RNA)
## get min, max and mean
mean_x <- round(mean(x = plot_data_df$nFeature_RNA), digits = 0)
mean_x
min_x <- round(min(x = plot_data_df$nFeature_RNA), digits = 0)
max_x <- round(max(x = plot_data_df$nFeature_RNA), digits = 0)

# make plot ---------------------------------------------------------------
## make color for the bar
fill_bar <- RColorBrewer::brewer.pal(n = 8, name = "Set1")[4]
p <- ggplot()
p <- p + geom_histogram(data = plot_data_df, mapping = aes(x = nFeature_RNA), color = "black", fill = fill_bar, binwidth = 200)
p <- p + ggtitle("Genes expressed per cell")
p <- p + xlab("Number of genes expressed") + ylab("Number of cells")
p <- p + scale_y_continuous(breaks = seq(0, 15000, 2500), labels = seq(0, 15000, 2500))
p <- p + xlim(c(0, 5000))
p <- p + annotate("text", label = paste0("Min = ", min_x, "\n",
                                         "Max = ", max_x, "\n",
                                         "Mean = ", mean_x), x = 1000, y = 12500)

p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(plot.title = element_text(hjust = 0.5, face = "bold"))
p

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "number_genes_per_cell_postQC.", "pdf")
pdf(file2write, width = 6, height = 4)
print(p)
dev.off()

