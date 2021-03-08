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
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups/20200917.v2/31Aliquot.Barcode2CellType.20200917.v2.tsv", data.table = F)
## input UMAP info per barcode
umap_df <- fread(input = "./Resources/Analysis_Results/data_summary/fetch_data/fetch_data_by_individual_sample//20200406.v1/individual_cluster_meta_data.20200406.v1.tsv", data.table = F)
## input id meta data table
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")

# plot for cell group----------------------------------------------------------
for (id_aliquot_tmp in c("CPT0001260013")) {
# for (id_aliquot_tmp in unique(umap_df$aliquot)) {
  aliquot_show <- idmetadata_df$Aliquot.snRNA.WU[idmetadata_df$Aliquot.snRNA == id_aliquot_tmp]
  
  plotdata_df <- umap_df %>%
    filter(aliquot == id_aliquot_tmp)
  plotdata_df <- merge(plotdata_df,
                       barcode2celltype_df %>%
                         select(individual_barcode, Cell_group15),
                       by.x = c("individual_barcode"), by.y = c("individual_barcode"), all.x = T)
  table(plotdata_df$Cell_group15)
  
  count_bycelltype_df <- plotdata_df %>%
    select(Cell_group15) %>%
    table() %>%
    as.data.frame() %>%
    rename(Cell_group15 = '.') %>%
    mutate(Keep = (Freq >= 50))
  plotdata_df <- plotdata_df %>%
    mutate(Cell_group = Cell_group15) %>%
    filter(Cell_group %in% count_bycelltype_df$Cell_group15[count_bycelltype_df$Keep])
  table(plotdata_df$Cell_group)
  
  ## make colors
  colors_nontumor <- colorblind_pal()((length(unique(plotdata_df$Cell_group)) - 2))
  # colors_nontumor <- Polychrome::dark.colors(n = (length(unique(plotdata_df$Cell_group)) - 1))[-2]
  colors_tumorcells <- cellgroup_colors[c("Tumor cells", "Unknown")]
  colors_cellgroup_tmp <- c(colors_tumorcells, colors_nontumor)
  names(colors_cellgroup_tmp) <- c("Tumor cells", "Unknown", unique(plotdata_df$Cell_group[!(plotdata_df$Cell_group %in% c("Tumor cells", "Unknown"))]))
  colors_cellgroup_tmp
  
  p <- ggplot()
  p <- p + geom_point(data = plotdata_df, 
                      mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_group),
                      alpha = 1, size = 0.2)
  p <- p + scale_color_manual(values = colors_cellgroup_tmp)
  p <- p + guides(colour = guide_legend(override.aes = list(size=3, fontsize = 20), nrow = ceiling(length(colors_cellgroup_tmp)/4), byrow = T))
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
  p <- p + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank())
  p <- p + theme(axis.text.y=element_blank(),
                 axis.ticks.y=element_blank())
  p <- p + ggtitle(label = paste0(aliquot_show, " Cell Types"))
  p <- p + theme(legend.position = "bottom")
  p
  file2write <- paste0(dir_out, aliquot_show, ".png")
  png(filename = file2write, width = 800, height = 1100, res = 150)
  print(p)
  dev.off()
  
  p <- ggplot()
  p <- p + geom_point(data = plotdata_df, 
                      mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_group),
                      alpha = 1, size = 0.2)
  p <- p + scale_color_manual(values = colors_cellgroup_tmp)
  p <- p + guides(colour = guide_legend(override.aes = list(size=5), nrow = ceiling(length(colors_cellgroup_tmp)/4), byrow = T))
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
  p <- p + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank())
  p <- p + theme(axis.text.y=element_blank(),
                 axis.ticks.y=element_blank())
  p <- p + theme(axis.title = element_text(size = 20))
  # p <- p + ggtitle(label = paste0(aliquot_show, " Cell Types"))
  p <- p + theme(legend.position = "bottom", legend.text = element_text(size = 20), legend.title = element_text(size = 25))
  p
  file2write <- paste0(dir_out, aliquot_show, ".pdf")
  pdf(file2write, width = 8, height = 9, useDingbats = F)
  print(p)
  dev.off()
}
