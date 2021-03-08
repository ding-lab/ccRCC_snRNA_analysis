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
version_tmp <- 1
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
table(plot_data_df$Cell_group) %>% length()
colors_tmp <- colors_cellgroup13
names(colors_tmp) <- names(table(plot_data_df$Cell_group))
p <- ggplot()
p <- p + geom_point_rast(data = plot_data_df, 
                         mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_group),
                         alpha = 1, size = 1, shape = 16)
p <- p + scale_color_manual(values = colors_tmp)
p <- p + theme_void()
p <- p + theme(legend.position="none")
p
## save as pdf
file2write <- paste0(dir_out, "rasterized.cellgroup_on_umap.", ".pdf")
pdf(file = file2write, width = 8, height = 8, useDingbats = F)
print(p)
dev.off()


# plot legend -------------------------------------------------------------
p <- ggplot()
p <- p + geom_point_rast(data = plot_data_df, 
                         mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_group),
                         alpha = 1, size = 1, shape = 16)
p <- p + scale_color_manual(values = colors_tmp)
p <- p + guides(colour = guide_legend(override.aes = list(size=3), nrow = 4))
p <- p + theme_void()
p <- p + theme(legend.position="bottom")
p
## save as pdf
file2write <- paste0(dir_out, "legend", ".pdf")
pdf(file = file2write, width = 9, height = 10, useDingbats = F)
print(p)
dev.off()

