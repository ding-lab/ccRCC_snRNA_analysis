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
# barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_corrected_by_individual_sample_inspection/20200828.v1/31Aliquot.Barcode2CellType.20200828.v1.tsv", data.table = F)
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_corrected_by_individual_sample_inspection/20201002.v1/31Aliquot.Barcode2CellType.20201002.v1.tsv", data.table = F)
## input UMAP info per barcode
umap_df <- fread(input = "./Resources/Analysis_Results/data_summary/fetch_data/fetch_data_by_individual_sample//20200406.v1/individual_cluster_meta_data.20200406.v1.tsv", data.table = F)
## input id meta data table
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")

# plot for cell group----------------------------------------------------------
for (id_aliquot_tmp in c("CPT0001260013")) {
# for (id_aliquot_tmp in c("CPT0075130004", "CPT0075140002", "CPT0075120002")) {
# for (id_aliquot_tmp in unique(umap_df$aliquot)) {
  aliquot_show <- idmetadata_df$Aliquot.snRNA.WU[idmetadata_df$Aliquot.snRNA == id_aliquot_tmp]
  
  plotdata_df <- umap_df %>%
    filter(aliquot == id_aliquot_tmp)
  plotdata_df <- merge(plotdata_df,
                       barcode2celltype_df %>%
                         filter(orig.ident == id_aliquot_tmp) %>%
                         select(individual_barcode, Cell_group.detailed),
                       by.x = c("individual_barcode"), by.y = c("individual_barcode"), all.x = T)
  plotdata_df <- plotdata_df %>%
    rename(Cell_group = Cell_group.detailed)
  p <- ggplot()
  p <- p + geom_point(data = plotdata_df, 
                      mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_group),
                      alpha = 1, size = 0.05)
  p <- p + scale_color_manual(values = cellgroup_colors[unique(plotdata_df$Cell_group)])
  p <- p + guides(colour = guide_legend(override.aes = list(size=3), ncol = ceiling(x = length(unique(plotdata_df$Cell_group))/2), byrow = T))
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
  p <- p + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank())
  p <- p + theme(axis.text.y=element_blank(),
                 axis.ticks.y=element_blank())
  p <- p + ggtitle(label = paste0(aliquot_show, " Cell Type"))
  p <- p + theme(legend.position = "top")
  p
  file2write <- paste0(dir_out, aliquot_show, ".png")
  png(filename = file2write, width = 800, height = 1000, res = 150)
  print(p)
  dev.off()
  
  p <- ggplot()
  p <- p + geom_point(data = plotdata_df, 
                      mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_group),
                      alpha = 1, size = 0.5)
  p <- p + scale_color_manual(values = cellgroup_colors[unique(plotdata_df$Cell_group)])
  p <- p + guides(colour = guide_legend(override.aes = list(size=5), nrow = ceiling(length(colors_cellgroup_tmp)/4), byrow = T))
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
  p <- p + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank())
  p <- p + theme(axis.text.y=element_blank(),
                 axis.ticks.y=element_blank())
  p <- p + theme(axis.title = element_text(size = 20))
  p <- p + theme(legend.position = "bottom", legend.text = element_text(size = 20), legend.title = element_text(size = 25))
  p
  file2write <- paste0(dir_out, aliquot_show, ".pdf")
  pdf(file2write, width = 8, height = 9, useDingbats = F)
  print(p)
  dev.off()
}
