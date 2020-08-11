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
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_to_all_cells/20200720.v1/30AliquotIntegration.Barcode2CellType.TumorManualCluster.20200720.v1.tsv", data.table = F)
## input UMAP info per barcode
umap_df <- fread(input = "./Resources/Analysis_Results/data_summary/fetch_data/fetch_data_by_individual_sample//20200406.v1/individual_cluster_meta_data.20200406.v1.tsv", data.table = F)

# plot for cell group----------------------------------------------------------
for (id_aliquot_tmp in unique(umap_df$aliquot)) {
  plotdata_df <- umap_df %>%
    filter(aliquot == id_aliquot_tmp)
  plotdata_df <- merge(plotdata_df,
                       barcode2celltype_df %>%
                         select(individual_barcode, Cell_group),
                       by.x = c("individual_barcode"), by.y = c("individual_barcode"), all.x = T)
  
  p <- ggplot()
  p <- p + geom_point(data = plotdata_df, 
                      mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_group),
                      alpha = 1, size = 0.05)
  p <- p + scale_color_manual(values = cellgroup_colors)
  p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
  p <- p + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank())
  p <- p + theme(axis.text.y=element_blank(),
                 axis.ticks.y=element_blank())
  p
  file2write <- paste0(dir_out, id_aliquot_tmp, ".cellgroup_on_umap.", run_id, ".png")
  png(filename = file2write, width = 1100, height = 800, res = 150)
  print(p)
  dev.off()
}
