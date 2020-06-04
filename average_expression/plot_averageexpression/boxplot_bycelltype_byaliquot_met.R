# Yige Wu @WashU Jun 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(ggridges)
library(viridis)
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input denpendencies -----------------------------------------------------
## input average expression by cell type by aliquot
avgexp_bycelltype_byaliquot_df <- fread(input = "./Resources/Analysis_Results/average_expression/averageexpression_bycelltypeshorter_byaliquot_on_katmai/20200411.v1/averageexpression_bycelltypeshorter.30_aliquot_integration.20200411.v1.tsv", data.table = F)
avgexp_bycelltype_byaliquot_df <- fread(input = "./Resources/Analysis_Results/average_expression/a", data.table = F)

## input id meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200413.v1/meta_data.20200413.v1.tsv", data.table = F)

# set plot gene list ------------------------------------------------------
gene_plot <- c("VHL")
colorpalette <- "Oranges"

# make matrix for heatmap body --------------------------------------------
plot_data_df <- avgexp_bycelltype_byaliquot_df %>%
  filter(V1 == gene_plot)
## melt
plot_data_long_df <- melt(plot_data_df)
## add cell type and aliquot info
plot_data_long_df <- plot_data_long_df %>%
  mutate(idaliquot_celltype = gsub(x = variable, pattern = "RNA.", replacement = "")) %>%
  mutate(celltype = str_split_fixed(string = idaliquot_celltype, pattern = "_", n = 2)[,2]) %>%
  mutate(id_aliquot = str_split_fixed(string = idaliquot_celltype, pattern = "_", n = 2)[,1])
plot_data_long_df$Sample_Type <- mapvalues(x = plot_data_long_df$id_aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Sample_Type))
plot_data_long_df <- plot_data_long_df %>%
  filter(Sample_Type != "Normal") %>%
  filter(!(celltype %in% c("Unknown", "Mixed.myeloid.lymphoid", "Macrophages.M2b", "CD4.CD8.proliferating", "TRM", "Myofibroblasts", "pDC", "Normal.epithelial.cells", "Basophils", "CD8..T.cells.exhausted")))
plot_data_long_df$id_aliquot_wu <- mapvalues(x = plot_data_long_df$id_aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
## cap value
cap_value <- quantile(x = plot_data_long_df$value, probs = 0.99)
cap_value
cap_value <- 0.6
plot_data_long_df <- plot_data_long_df %>%
  mutate(value_capped = ifelse(value > cap_value, NA, value))
color_fill <- RColorBrewer::brewer.pal(n = 1, name = colorpalette)[3]
# plot ridge plot ---------------------------------------------------------
p <- ggplot()
p <- p + geom_boxplot(data = plot_data_long_df, mapping = aes(x = value_capped, y = celltype), alpha = 0.8, fill = color_fill)
p <- p + theme_bw()
p <- p + theme(panel.spacing = unit(0.1, "lines"),
               strip.text.x = element_text(size = 10, colour = "black"),
               axis.text.y = element_text(size = 12))
p <- p + xlab("Average Expression Value by Cell Type")
p <- p + ylab("Cell Type")
p <- p + labs(fill = "Average expression value\nby cell type (Normalized\nCapped at 99% quantile)")
p

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, gene_plot, ".expression.", "by_celltype.", run_id, ".png")
png(filename = file2write, width = 1500, height = 1200, res = 150)
print(p)
dev.off()

file2write <- paste0(dir_out, gene_plot, ".expression.", "by_celltype.", run_id, ".pdf")
pdf(file = file2write, width = 8, height = 6)
print(p)
dev.off()

