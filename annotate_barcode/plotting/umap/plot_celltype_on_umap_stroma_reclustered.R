# Yige Wu @WashU Apr 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
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
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_corrected_by_individual_sample_inspection/20201002.v1/31Aliquot.Barcode2CellType.20201002.v1.tsv", data.table = F)
## input UMAP info per barcode
umap_df <- fread(input = "./Resources/Analysis_Results/recluster/recluster_cell_groups_in_integrated_data/fetch_data/fetch_data_stroma_reclustered/20201005.v1/stroma_reclustered.metadata.20201005.v1.tsv", data.table = F)

# make data frame for plotting---------------------------------------------------------
plotdata_df <- umap_df %>%
  select(orig.ident, original_barcode, UMAP_1, UMAP_2)
plotdata_df <- merge(plotdata_df,
                     barcode2celltype_df %>%
                       mutate(Cell_type = Cell_type.shorter) %>%
                       select(orig.ident, individual_barcode, Cell_type),
                     by.x = c("orig.ident", "original_barcode"), by.y = c("orig.ident", "individual_barcode"), all.x = T)


# plot --------------------------------------------------------------------
unique(plotdata_df$Cell_type)
p <- ggplot()
p <- p + geom_point(data = plotdata_df, 
                    mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_type),
                    alpha = 1, size = 0.05)
p <- p + scale_color_manual(values = stroma_colors)
p <- p + guides(colour = guide_legend(override.aes = list(size=3)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, aliquot_show, ".png")
png(filename = file2write, width = 1200, height = 1000, res = 150)
print(p)
dev.off()
