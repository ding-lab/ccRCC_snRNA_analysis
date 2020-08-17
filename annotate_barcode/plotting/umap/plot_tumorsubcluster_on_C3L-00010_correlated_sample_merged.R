# Yige Wu @WashU Aug 2020

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
## input umap data
umap_df <- fread(data.table = F, input = "./Resources/Analysis_Results/fetch_data/fetchdata_C3L-00010_C3L-00416_C3L-00583_C3N-00242_merged_on_katmai/20200817.v1/C3L-00010_C3L-00416_C3L-00583_C3N-00242.Metadata.20200817.v1.tsv")
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")

# make plot data ----------------------------------------------------------
plotdata_df <- umap_df %>%
  mutate(individual_barcode = str_split_fixed(string = barcode, pattern = "_", n = 2)[,1])
## map readble id
plotdata_df$Id_Aliquot_WU <- mapvalues(x = plotdata_df$orig.ident, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
unique(plotdata_df$Id_Aliquot_WU)
## map cell type and tumor subcluster id
plotdata_df <- merge(plotdata_df, barcode2celltype_df, by = c("orig.ident", "individual_barcode"), all.x = T)
plotdata_df <- plotdata_df %>%
  mutate(Name_Cluster = paste0("C", (Id_TumorManualCluster + 1)))

# plot by csample--------------------------------------------------------------------
p <- ggplot()
p <- p + geom_point(data = plotdata_df, 
                    mapping = aes(x = UMAP_1, y = UMAP_2, color = Id_Aliquot_WU),
                    alpha = 1, size = 0.05)
p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p

# plot by cell type--------------------------------------------------------------------
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

# plot by tumor subcluster id--------------------------------------------------------------------
for (id_sample in unique(plotdata_df$Id_Aliquot_WU)) {
  plotdata_sub_df <- plotdata_df %>%
    filter(Id_Aliquot_WU == id_sample)
  p <- ggplot()
  p <- p + geom_point(data = plotdata_sub_df, 
                      mapping = aes(x = UMAP_1, y = UMAP_2, color = Name_Cluster),
                      alpha = 1, size = 0.05)
  p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
  p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                 panel.background = element_blank(), axis.line = element_line(colour = "black"))
  p <- p + theme(axis.text.x=element_blank(),
                 axis.ticks.x=element_blank())
  p <- p + theme(axis.text.y=element_blank(),
                 axis.ticks.y=element_blank())
  p
  ## save output
  file2write <- paste0(dir_out, id_sample, ".png")
  png(file2write, width = 1000, height = 800, res = 150)
  print(p)
  dev.off()
}


