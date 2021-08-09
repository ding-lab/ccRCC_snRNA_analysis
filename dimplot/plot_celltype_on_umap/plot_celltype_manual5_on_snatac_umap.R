# Yige Wu @WashU May 2020
## plot cell type on integration UMAP

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
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
names(colors_cellgroup14)
colors_tmp1 <- colors_cellgroup14[c("Tumor cells", "CD4+ T-cells", "CD8+ T-cells", "Macrophages", "NK cells", "DC", "Fibroblasts")]
colors_tmp2 <- c(colors_cellgroup14[c("Transitional cells", "Normal epithelial cells", "Myofibroblasts", "Immune others", "B-cells")], 
                 Polychrome::palette36.colors(n = 36)[c("Vivid_Violet","Light_Olive_Brown", "Very_Light_Blue")], "grey80")
names(colors_tmp2) <- c("EMT tumor cells", "PT", "Loop of Henle", "Distal convoluted tubule", 'Principle cells', 
                        "Intercalated cells", "Podocytes", "Endothelial cells", "Unknown")
colors_tmp <- c(colors_tmp1, colors_tmp2)
swatch(colors_tmp)

p <- ggplot()
p <- p + geom_point_rast(data = plot_data_df, 
                         mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_group),
                         alpha = 1, size = 0.05, shape = 16)
p <- p + scale_color_manual(values = colors_tmp)
p <- p + guides(colour = guide_legend(override.aes = list(size=5), nrow = 4))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.line=element_blank(),axis.text.x=element_blank(),
               axis.text.y=element_blank(),axis.ticks=element_blank(),
               axis.title.x=element_blank(),
               axis.title.y=element_blank(), 
               panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
               panel.grid.minor=element_blank(),plot.background=element_blank())
p <- p + theme(legend.position="bottom")
p
## save as pdf
file2write <- paste0(dir_out, "cellgroup_on_umap.", ".pdf")
pdf(file = file2write, width = 8, height = 10, useDingbats = F)
print(p)
dev.off()
## save as png
file2write <- paste0(dir_out, "cellgroup_on_umap.", ".png")
png(filename = file2write, width = 1000, height = 1100, res = 150)
print(p)
dev.off()

#