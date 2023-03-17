# Yige Wu @WashU May 2020
## plot cell type on integration UMAP

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
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
## make output directory
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# input  ------------------------------------------------------
## input expression data
exp_df <- readRDS(file = "./Resources/Analysis_Results/fetch_data/extract_exp/extract_selected_markers__sct_data_bybarcode/20230130.v1/tumor_specific_marker_expression_bybarcode.RDS")
## input cell type per barcode table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups_35aliquots/20210802.v1/35Aliquot.Barcode2CellType.20210802.v1.tsv", data.table = F)
## input UMAP info per barcode
integrated_umap_df <- fread(input = "./Resources/Analysis_Results/fetch_data/fetchdata_34_ccRCC_samples_merged_katmai/20211005.v1/ccRCC.34Sample.Merged.Metadata.20211005.v1.tsv", data.table = F)

# preprocess cell type annotation for each barcode ----------------------------------------------------------
barcode_df <- merge(integrated_umap_df %>%
                        mutate(individual_barcode = str_split_fixed(string = barcode, pattern = "_", n = 2)[,1]),
                      barcode2celltype_df %>%
                        mutate(Cell_group = Cell_group_w_epithelialcelltypes) %>%
                        select(orig.ident, individual_barcode, Cell_group),
                      by = c("orig.ident", "individual_barcode"), all.x = T
)
barcode_df <- barcode_df %>%
  arrange(desc(Cell_group))
barcode_df <- rbind(barcode_df[barcode_df$Cell_group %in% c("Unknown", "Immune others"),], barcode_df[!(barcode_df$Cell_group %in% c("Unknown", "Immune others")),])
table(barcode_df$Cell_group)

# specify genes to plot ---------------------------------------------------
gene_x <- "EPAS1"
gene_y <- "CD70"
# gene_y <- "EPAS1"

# make plot data ----------------------------------------------------------
exp_df$barcode <- rownames(exp_df)
exp_df$value_x <- exp_df[,gene_x]
exp_df$value_y <- exp_df[,gene_y]
plot_data_df <- merge(x = exp_df[, c("barcode", "value_x", "value_y")],
                      y = barcode_df[, c("barcode", "Cell_group")], by.x = c("barcode"), by.y = c("barcode"), all.x = T)
plot_filtereddata_df <- plot_data_df %>%
  filter(value_x > 0 & value_y > 0)

# make plot - tumor cells with expression--------------------------------------------------------------
plot_tmpdata_df <- plot_filtereddata_df %>%
  filter(Cell_group == "Tumor cells")
summary(plot_tmpdata_df$value_x)
table(plot_tmpdata_df$value_x)
table(plot_tmpdata_df$value_y)

plot_tmpdata_df <- plot_tmpdata_df %>%
  mutate(value_x_bin = cut(value_x, breaks = c(0, 0.70, 3.3)))
table(plot_tmpdata_df %>%
        select(value_x, value_x_bin))
p <- ggplot(data = plot_tmpdata_df, mapping = aes(x = value_x_bin, y = value_y))
p <- p + geom_violin(alpha = 0.5)
p <- p + geom_boxplot(alpha = 0.5)
# p <- p + geom_jitter_rast(alpha = 0.8, shape = 16)
# Add correlation coefficient
p <- p + stat_cor(method = "pearson", label.x = 1, label.y = 1, size = 6)
p <- p + xlab(paste0(gene_x, " expression"))
p <- p + ylab(paste0(gene_y, " expression"))
p <- p + theme_classic(base_size = 16)
p <- p + theme(axis.text = element_text(color = "black", size = 16),
               axis.title = element_text(color = "black", size = 16))
p
## save as png
file2write <- paste0(dir_out, gene_x, ".vs.", gene_y, ".", "Tumorcells", "_w_expression", ".png")
png(filename = file2write, width = 1000, height = 1100, res = 150)
print(p)
dev.off()

# make plot - tumor cells with expression--------------------------------------------------------------
plot_tmpdata_df <- plot_data_df %>%
  filter(Cell_group == "Tumor cells")
summary(plot_tmpdata_df$value_x)
table(plot_tmpdata_df$value_x)
table(plot_tmpdata_df$value_y)

plot_tmpdata_df <- plot_tmpdata_df %>%
  mutate(value_x_bin = cut(value_x, breaks = c(-1, 0.1, 0.69, 0.70, 3.3)))
table(plot_tmpdata_df %>%
        select(value_x, value_x_bin))
p <- ggplot(data = plot_tmpdata_df, mapping = aes(x = value_x_bin, y = value_y))
p <- p + geom_violin(alpha = 0.5)
p <- p + geom_boxplot(alpha = 0.5)
# p <- p + geom_jitter_rast(alpha = 0.8, shape = 16)
# Add correlation coefficient
p <- p + stat_cor(method = "pearson", label.x = 1, label.y = 1, size = 6)
p <- p + xlab(paste0(gene_x, " expression"))
p <- p + ylab(paste0(gene_y, " expression"))
p <- p + theme_classic(base_size = 16)
p <- p + theme(axis.text = element_text(color = "black", size = 16),
               axis.title = element_text(color = "black", size = 16))
p
## save as png
file2write <- paste0(dir_out, gene_x, ".vs.", gene_y, ".", "Tumorcells", ".png")
png(filename = file2write, width = 1000, height = 1100, res = 150)
print(p)
dev.off()

