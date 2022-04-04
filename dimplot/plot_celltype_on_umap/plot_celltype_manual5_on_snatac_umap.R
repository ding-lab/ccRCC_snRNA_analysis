# Yige Wu @WashU May 2020
## plot cell type on integration UMAP

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(ggrastr)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input cell type per barcode table
# barcode2celltype_df <- fread(input = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/Barcode_Annotation/UMAP/UMAP_data_13_snATAC_Normal_epithelialCells_reclustered.20201209.tsv", data.table = F)
barcode2celltype_df <- fread(input = "./Resources/snATAC_Processed_Data/Barcode_Annotation/28_snATAC.ManualReviwed.UMAP_data.20210709.tsv", data.table = F)

# plot for cell group----------------------------------------------------------
table(barcode2celltype_df$cell_type)
plot_data_df <- barcode2celltype_df %>%
  mutate(Cell_group = ifelse(cell_type %in% c("EMT tumor cells", "Tumor"), "Tumor cells", cell_type))
  # mutate(Cell_group = ifelse(cell_type_manual_5 == "EMT tumor cells", "Tumor cells", cell_type_manual_5))
table(plot_data_df$Cell_group)
## make colors
colors_tmp1 <- colors_cellgroup14[c("Tumor cells", "CD4+ T-cells", "CD8+ T-cells", "Macrophages", "NK cells", "DC", "Fibroblasts", "Myofibroblasts",  "B-cells")]
colors_tmp2 <- c(colors_cellgroup14[c("Normal epithelial cells", "EMT tumor cells", "Immune others")],
                 Polychrome::palette36.colors(n = 36)[c("Vivid_Yellow_Green", "Vivid_Violet","Light_Olive_Brown", "Very_Light_Blue")], 
                 "grey80", "grey50")
names(colors_tmp2) <- c("PT", "Loop of Henle", "Distal convoluted tubule", 
                        'Principle cells', "Intercalated cells", "Podocytes", "Endothelial cells", 
                        "Unknown", "Immune others")
colors_cellgroup <- c(colors_tmp1, colors_tmp2)
# swatch(colors_cellgroup)
# swatch(Polychrome::palette36.colors(n = 36))

p <- ggplot()
p <- p + geom_point_rast(data = plot_data_df, 
                         mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_group),
                         alpha = 1, size = 0.05, shape = 16)
p <- p + scale_color_manual(values = colors_cellgroup)
p <- p + guides(colour = guide_legend(override.aes = list(size=4)))
p <- p + theme_void()
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p <- p + theme(legend.position="bottom", aspect.ratio=1)
p
## save as pdf
file2write <- paste0(dir_out, "cellgroup_on_umap.", ".pdf")
pdf(file = file2write, width = 8, height = 9, useDingbats = F)
print(p)
dev.off()
## save as png
file2write <- paste0(dir_out, "cellgroup_on_umap.", ".png")
png(filename = file2write, width = 1000, height = 1100, res = 150)
print(p)
dev.off()

#