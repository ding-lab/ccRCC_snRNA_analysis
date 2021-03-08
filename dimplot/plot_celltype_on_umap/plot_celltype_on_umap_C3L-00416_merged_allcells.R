# Yige Wu @WashU Apr 2020

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
## input id meta data table
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")

# input seurat object, and get umap info -----------------------------------------------------
aliquot_show <- "C3L-00416 Merged"
srat <- readRDS(file = "./Data Freezes/V1/snRNA/Merged_Seurat_Objects/C3L-00416.Tumor_Segments.Merged.20200319.v1.RDS")
umap_df <- FetchData(object = srat, vars = c("orig.ident", "UMAP_1", "UMAP_2"))
umap_df$barcode_integrated <- rownames(umap_df)

# make colors for the cell types ------------------------------------------
colors_plot <- c(colors_cellgroup14, Polychrome::palette36.colors(n = 36)[c(13,18)])
names(colors_plot) <- c(names(colors_cellgroup14), "Mixed myeloid/lymphoid", "CD4/CD8 proliferating")

# plot with legend----------------------------------------------------------
plotdata_df <- umap_df %>%
  mutate(individual_barcode = str_split_fixed(string = barcode_integrated, pattern = "_", n = 3)[,1])
plotdata_df <- merge(plotdata_df,
                     barcode2celltype_df %>%
                       mutate(Cell_group = ifelse(Cell_type.detailed %in% c("Mixed myeloid/lymphoid", "CD4/CD8 proliferating"), Cell_type.detailed, Cell_group14_w_transitional)) %>%
                       select(orig.ident, individual_barcode, Cell_group, Cell_type.detailed),
                     by.x = c("orig.ident", "individual_barcode"), by.y = c("orig.ident", "individual_barcode"), all.x = T)
p <- ggplot()
p <- p + geom_point_rast(data = plotdata_df, 
                    mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_group),
                    alpha = 1, size = 0.05)
p <- p + scale_color_manual(values = colors_plot)
p <- p + guides(colour = guide_legend(override.aes = list(size=3)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p
file2write <- paste0(dir_out, aliquot_show, ".png")
png(filename = file2write, width = 1200, height = 1000, res = 150)
print(p)
dev.off()
file2write <- paste0(dir_out, aliquot_show, ".withlegend.pdf")
pdf(file2write, width = 12, height = 10, useDingbats = F)
print(p)
dev.off()

# plot without legend----------------------------------------------------------
plotdata_df <- umap_df %>%
  mutate(individual_barcode = str_split_fixed(string = barcode_integrated, pattern = "_", n = 3)[,1])
plotdata_df <- merge(plotdata_df,
                     barcode2celltype_df %>%
                       mutate(Cell_group = ifelse(Cell_type.detailed %in% c("Mixed myeloid/lymphoid", "CD4/CD8 proliferating"), Cell_type.detailed, Cell_group14_w_transitional)) %>%
                       select(orig.ident, individual_barcode, Cell_group, Cell_type.detailed),
                     by.x = c("orig.ident", "individual_barcode"), by.y = c("orig.ident", "individual_barcode"), all.x = T)
p <- ggplot()
p <- p + geom_point_rast(data = plotdata_df, 
                         mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_group),
                         alpha = 1, size = 0.05)
p <- p + scale_color_manual(values = colors_plot)
p <- p + guides(colour = guide_legend(override.aes = list(size=3)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p <- p + theme(legend.position = "none")
p
file2write <- paste0(dir_out, aliquot_show, ".nolegend.pdf")
pdf(file2write, width = 4, height = 4, useDingbats = F)
print(p)
dev.off()

