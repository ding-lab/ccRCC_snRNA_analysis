# Yige Wu @WashU May 2020
## plot cell type on integration UMAP

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
if (!("ggrastr" %in% installed.packages()[,1])) {
  install.packages('ggrastr')
}
library(ggrastr)
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

# input dependencies ------------------------------------------------------
## input cell type per barcode table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups/20201002.v1/31Aliquot.Barcode2CellType.20201002.v1.tsv", data.table = F)
## input UMAP info per barcode
integrated_umap_df <- fread(input = "./Resources/Analysis_Results/integration/30_aliquot_integration/fetch_data/20200212.v3/30_aliquot_integration.20200212.v3.umap_data.tsv", data.table = F)

# plot for cell group----------------------------------------------------------
plot_data_df <- merge(integrated_umap_df,
                      barcode2celltype_df %>%
                        filter(!is.na(integrated_barcode)) %>%
                        mutate(Cell_group = Cell_group7) %>%
                        select(integrated_barcode, Cell_group),
                      by.x = c("barcode"), by.y = c("integrated_barcode"), all.x = T)
plot_data_df$Cell_group[plot_data_df$Cell_group == "Tumor-like cells"] <- "Tumor cells"
plot_data_df$Cell_group[plot_data_df$Cell_group == "Transitional cells"] <- "Tumor cells"
plot_data_df$Cell_group[plot_data_df$Cell_group == "Normal-like cells"] <- "Normal epithelial cells"
unique(plot_data_df$Cell_group)
p <- ggplot()
p <- p + geom_point_rast(data = plot_data_df, 
                    mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_group),
                    alpha = 1, size = 0.05, shape = 16)
p <- p + scale_color_manual(values = cellgroup_colors)
p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p <- p + theme(legend.position="bottom")
p
## save as pdf
file2write <- paste0(dir_out, "cellgroup_on_umap.", run_id, ".pdf")
pdf(file = file2write, width = 8, height = 9, useDingbats = F)
print(p)
dev.off()
## save as png
file2write <- paste0(dir_out, "cellgroup_on_umap.", run_id, ".png")
png(filename = file2write, width = 1000, height = 1100, res = 150)
print(p)
dev.off()