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
umap_df <- fread(data.table = F, input = "./Resources/Analysis_Results/fetch_data/fetchdata_tumorcells_C3L-00010_C3L-00416_C3L-00583_C3N-00242_merged_on_katmai/20200817.v1/C3L-00010_C3L-00416_C3L-00583_C3N-00242.Metadata.20200817.v1.tsv")
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

# plot by sample--------------------------------------------------------------------
## make color palette
colors_sample <- RColorBrewer::brewer.pal(n = 5, name = "Set1")
names(colors_sample) <- unique(plotdata_df$Id_Aliquot_WU)
p <- ggplot()
p <- p + geom_point(data = plotdata_df, 
                    mapping = aes(x = UMAP_1, y = UMAP_2, color = Id_Aliquot_WU),
                    alpha = 0.8, size = 0.05)
p <- p + scale_color_manual(values = colors_sample)
p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p <- p + theme(legend.position = "top")
p
## save output
file2write <- paste0(dir_out, "bysample", ".png")
png(file2write, width = 1000, height = 1000, res = 150)
print(p)
dev.off()
file2write <- paste0(dir_out, "bysample", ".pdf")
pdf(file2write, width = 6, height = 7, useDingbats = F)
print(p)
dev.off()
file2write <- paste0(dir_out, "bysample.legend", ".pdf")
pdf(file2write, width = 9, height = 10, useDingbats = F)
print(p)
dev.off()

# plot by cluster--------------------------------------------------------------------
## make factor
plotdata_df$Reclustered_Cluster_Number <- factor(plotdata_df$ident)
## make color
colors_reclustered_cluster <- c(RColorBrewer::brewer.pal(n = 8, name = "Dark2"), RColorBrewer::brewer.pal(n = 5, name = "Set1"))
names(colors_reclustered_cluster) <- 0:12
p <- ggplot()
p <- p + geom_point(data = plotdata_df, 
                    mapping = aes(x = UMAP_1, y = UMAP_2, color = Reclustered_Cluster_Number),
                    alpha = 1, size = 0.05)
p <- p + scale_color_manual(values = colors_reclustered_cluster)
p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p
## save output
file2write <- paste0(dir_out, "bycluster", ".png")
png(file2write, width = 1200, height = 800, res = 150)
print(p)
dev.off()

# plot by tumor subcluster id--------------------------------------------------------------------
## make color
## make color for the cluster name suffix
unique(names_cluster_suffix)
colors_clustername <- c(RColorBrewer::brewer.pal(n = 8, name = "Dark2")[1:7], RColorBrewer::brewer.pal(n = 3, name = "Set1")[1], "grey50")
names(colors_clustername) <- c(paste0("C", 1:8), "other")
swatch(colors_clustername)

for (id_sample in unique(plotdata_df$Id_Aliquot_WU)) {
  plotdata_sub_df <- plotdata_df %>%
    mutate(Name_Cluster_tmp = ifelse(Id_Aliquot_WU == id_sample, Name_Cluster, "other"))
  p <- ggplot()
  p <- p + geom_point(data = plotdata_sub_df, 
                      mapping = aes(x = UMAP_1, y = UMAP_2, color = Name_Cluster_tmp),
                      alpha = 0.8, size = 0.05)
  p <- p + scale_color_manual(values = colors_clustername)
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


