# Manuscript:   Epigenetic and transcriptomic characterization reveals progression markers and essential pathways in clear cell renal cell carcinoma
# Author:       Yige Wu
# Description:  UMAP visualization of merged snRNA-seq data
# Section:      Results - Overview of clinical features and datasets
#=======================================================================================

# Load libraries  and set working directory  -----------------------------------
## load libraries
packages = c(
  "data.table",
  "dplyr",
  "ggrastr",
  "ggplot2"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    install.packages(pkg_name_tmp, dependencies = T)
  }
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg_name_tmp)
  }
  library(package = pkg_name_tmp, character.only = T)
}
## set working directory to current file location
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# input -------------------------------------------------------------------
plot_data_df <- fread(data.table = F, input = "../../data/snRNA.UMAPCoordinate.ClusterID.CellType.ByBarcode.tsv.gz")

# process data ------------------------------------------------------------
plot_data_df <- plot_data_df %>%
  mutate(Cell_group = Cell_group_w_epithelialcelltypes) %>%
  arrange(desc(Cell_group))
plot_data_df <- rbind(plot_data_df[plot_data_df$Cell_group %in% c("Unknown", "Immune others"),], plot_data_df[!(plot_data_df$Cell_group %in% c("Unknown", "Immune others")),])

# make colors -------------------------------------------------------------
colors_cellgroup <- c("#E7298A", "#E69F00", "#56B4E9", "#F0E442", "#D55E00", "#0072B2", "#FB9A99", "#B2DF8A", "#000000",
                      "#1B9E77", "#B15928", "#7570B3", "#90AD1C", "#AA0DFE", "#85660D", "#BDCDFF", "grey80", "grey50")
names(colors_cellgroup) <- c("Tumor cells", "CD4+ T-cells", "CD8+ T-cells", "Macrophages", "NK cells", "DC", "Fibroblasts", "Myofibroblasts",  "B-cells",
                             "Proximal tubule", "Loop of Henle", "Distal convoluted tubule", 
                             'Principle cells', "Intercalated cells", "Podocytes", "Endothelial cells", 
                             "Unknown", "Immune others")

# plot --------------------------------------------------------------------
p <- ggplot()
p <- p + geom_point_rast(data = plot_data_df, 
                         mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_group),
                         alpha = 1, size = 0.1, shape = 16)
p <- p + scale_color_manual(values = colors_cellgroup)
p <- p + guides(colour = guide_legend(override.aes = list(size=4), title = NULL))
p <- p + theme_void()
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_blank())
# axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p <- p + theme(legend.position="bottom", aspect.ratio=1)

#Print extended figure 2b
dir_out <- paste0(dir_base, "outputs/"); dir.create(dir_out)
file2write <- paste0(dir_out,"F1b_UMAP_snRNA.pdf")
pdf(file = file2write, width = 8, height = 9, useDingbats = F)
print(p)
dev.off()

