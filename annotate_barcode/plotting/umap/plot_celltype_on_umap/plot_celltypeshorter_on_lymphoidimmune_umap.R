# Yige Wu @WashU Dec 2020

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
version_tmp <- 4
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input cell type per barcode table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups/20201130.v1/31Aliquot.Barcode2CellType.20201130.v1.tsv", data.table = F)
## input UMAP info per barcode
integrated_umap_df <- fread(input = "./Resources/Analysis_Results/fetch_data/fetchdata_lymphoid_immunecells_on_katmai/20201210.v1/Metadata.20201210.v1.tsv", data.table = F)

# plot for cell group----------------------------------------------------------
plot_data_df <- merge(integrated_umap_df %>%
                        mutate(individual_barcode = str_split_fixed(string = barcode, pattern = "_", n = 2)[,1]),
                      barcode2celltype_df %>%
                        mutate(Cell_group = Cell_type.shorter) %>%
                        select(orig.ident, individual_barcode, Cell_group),
                      by = c("orig.ident", "individual_barcode"), all.x = T
)
table(plot_data_df$Cell_group)
p <- ggplot()
p <- p + geom_point_rast(data = plot_data_df, 
                         mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_group),
                         alpha = 1, size = 1, shape = 16)
# p <- p + scale_color_manual(values = cellgroup_colors)
p <- p + guides(colour = guide_legend(override.aes = list(size=5), nrow = 1))
p <- p + theme_void()
# p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#                panel.background = element_blank(), axis.line = element_line(colour = "black"))
# p <- p + theme(axis.text.x=element_blank(),
#                axis.ticks.x=element_blank())
# p <- p + theme(axis.text.y=element_blank(),
#                axis.ticks.y=element_blank())
p <- p + theme(legend.position="none")
p
## save as pdf
file2write <- paste0(dir_out, "rasterized.cellgroup_on_umap.", ".pdf")
pdf(file = file2write, width = 8, height = 8, useDingbats = F)
print(p)
dev.off()
