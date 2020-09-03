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
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input cell type per barcode table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_corrected_by_individual_sample_inspection/20200828.v1/31Aliquot.Barcode2CellType.20200828.v1.tsv", data.table = F)
## input UMAP info per barcode
umap_df <- fread(input = "./Resources/Analysis_Results/data_summary/fetch_data/fetch_data_by_individual_sample/20200717.v1/Barcode2MetaData.20200717.v1.tsv", data.table = F)
## input id meta data table
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")

# plot for cell group----------------------------------------------------------
for (id_aliquot_tmp in unique(umap_df$aliquot)) {
  aliquot_show <- idmetadata_df$Aliquot.snRNA.WU[idmetadata_df$Aliquot.snRNA == id_aliquot_tmp]
  plotdata_df <- umap_df %>%
    filter(aliquot == id_aliquot_tmp)
  plotdata_df <- merge(plotdata_df,
                       barcode2celltype_df %>%
                         select(individual_barcode, Cell_group.detailed, Cell_type.shorter),
                       by.x = c("individual_barcode"), by.y = c("individual_barcode"), all.x = T)
  
  count_bycelltype_bycluster_df <- plotdata_df %>%
    select(Cell_type.shorter, seurat_cluster_id) %>%
    table() %>%
    data.frame() %>%
    mutate(Id_Cluster_CellType = paste0("C", seurat_cluster_id, "_", Cell_type.shorter)) %>%
    mutate(Keep = (Freq >= 30))
  plotdata_df <- plotdata_df %>%
    mutate(Id_Cluster_CellType = paste0("C", seurat_cluster_id, "_", Cell_type.shorter))
  plotdata_df$Keep <- mapvalues(x = plotdata_df$Id_Cluster_CellType, from = count_bycelltype_bycluster_df$Id_Cluster_CellType, to = as.vector(count_bycelltype_bycluster_df$Keep))
  plotdata_df <- plotdata_df %>%
    filter(Keep == T)
  
  p <- ggplot()
  # p <- p + geom_violin(data = plotdata_df, 
  #                     mapping = aes(x = Cell_group.detailed, y = nCount_RNA, fill = Cell_type.shorter),
  #                     alpha = 0.8, size = 0.05)
  p <- p + geom_boxplot(data = plotdata_df, 
                        mapping = aes(x = Cell_type.shorter, y = nCount_RNA, fill = Cell_group.detailed), width = 0.4)
  p <- p + scale_fill_manual(values = cellgroup_colors)
  p <- p + facet_grid(~seurat_cluster_id,
                      scales = "free_x", space = "free", shrink = T)
  p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
  p <- p + theme(axis.text.x=element_text(angle = 90))
  p <- p + ggtitle(label = paste0(aliquot_show, " Gene Count By Cluster By Major Cell Type"))
  # p <- p + theme(axis.text.y=element_blank(),
  #                axis.ticks.y=element_blank())
  p
  file2write <- paste0(dir_out, aliquot_show, ".png")
  png(filename = file2write, width = 1000, height = 600, res = 150)
  print(p)
  dev.off()
}
