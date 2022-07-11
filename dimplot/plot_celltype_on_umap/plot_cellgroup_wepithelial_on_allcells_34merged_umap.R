# Yige Wu @WashU May 2020
## plot cell type on integration UMAP

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(ggrastr)

# input dependencies ------------------------------------------------------
## input cell type per barcode table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups_35aliquots/20210802.v1/35Aliquot.Barcode2CellType.20210802.v1.tsv", data.table = F)
## input UMAP info per barcode
integrated_umap_df <- fread(input = "./Resources/Analysis_Results/fetch_data/fetchdata_34_ccRCC_samples_merged_katmai/20211005.v1/ccRCC.34Sample.Merged.Metadata.20211005.v1.tsv", data.table = F)

# make plot data----------------------------------------------------------
plot_data_df <- merge(integrated_umap_df %>%
                        mutate(individual_barcode = str_split_fixed(string = barcode, pattern = "_", n = 2)[,1]),
                      barcode2celltype_df %>%
                        mutate(Cell_group = Cell_group_w_epithelialcelltypes) %>%
                        select(orig.ident, individual_barcode, Cell_group),
                      by = c("orig.ident", "individual_barcode"), all.x = T
)
plot_data_df <- plot_data_df %>%
  arrange(desc(Cell_group))
plot_data_df <- rbind(plot_data_df[plot_data_df$Cell_group %in% c("Unknown", "Immune others"),], plot_data_df[!(plot_data_df$Cell_group %in% c("Unknown", "Immune others")),])
table(plot_data_df$Cell_group)

# make colors -------------------------------------------------------------
colors_tmp1 <- colors_cellgroup14[c("Tumor cells", "CD4+ T-cells", "CD8+ T-cells", "Macrophages", "NK cells", "DC", "Fibroblasts", "Myofibroblasts",  "B-cells")]
colors_tmp2 <- c(colors_cellgroup14[c("Normal epithelial cells", "EMT tumor cells", "Immune others")],
                 Polychrome::palette36.colors(n = 36)[c("Vivid_Yellow_Green", "Vivid_Violet","Light_Olive_Brown", "Very_Light_Blue")], 
                 "grey80", "grey50")
names(colors_tmp2) <- c("Proximal tubule", "Loop of Henle", "Distal convoluted tubule", 
                        'Principle cells', "Intercalated cells", "Podocytes", "Endothelial cells", 
                        "Unknown", "Immune others")
colors_cellgroup <- c(colors_tmp1, colors_tmp2)
# swatch(colors_cellgroup)
# swatch(Polychrome::palette36.colors(n = 36))

# make plots --------------------------------------------------------------
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
p


# save output -------------------------------------------------------------
source("./ccRCC_snRNA_analysis/functions.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
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


# plot for response -------------------------------------------------------
p <- ggplot()
p <- p + geom_point_rast(data = plot_data_df, 
                         mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_group),
                         alpha = 1, size = 0.1, shape = 16)
p <- p + scale_color_manual(values = colors_cellgroup)
p <- p + theme_void()
p <- p + guides(colour = guide_legend(override.aes = list(size=3), title = NULL, ncol = 3, label.theme = element_text(size = 14)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(),
               axis.line = element_blank())
# axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p <- p + theme(legend.position="bottom")
p
## save as pdf
file2write <- paste0(dir_out, "cellgroup_on_umap2.", ".pdf")
pdf(file = file2write, width = 6, height = 7, useDingbats = F)
print(p)
dev.off()