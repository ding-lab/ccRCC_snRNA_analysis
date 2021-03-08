# Yige Wu @WashU May 2020
## plot cell type on integration UMAP

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 3
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input cell type per barcode table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_to_all_cells/20200626.v1/30AliquotIntegration.Barcode2CellType.TumorManualCluster.20200626.v1.tsv", data.table = F)
## input UMAP info per barcode
integrated_umap_df <- fread(input = "./Resources/Analysis_Results/integration/30_aliquot_integration/fetch_data/20200212.v3/30_aliquot_integration.20200212.v3.umap_data.tsv", data.table = F)
## input bulk genomics profile
bulkprofile_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/other/merge_bulk_events/20200512.v1/merged_bulk_events.20200512.v1.tsv")
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200505.v1/meta_data.20200505.v1.tsv")

# plot tumor cells by SETD2----------------------------------------------------------
## map cell type to each barcode
plot_data_df <- merge(integrated_umap_df,
                      barcode2celltype_df %>%
                        select(integrated_barcode, Cell_group),
                      by.x = c("barcode"), by.y = c("integrated_barcode"), all.x = T
)
## filter to just tumor cells
plot_data_df <- plot_data_df %>%
  filter(Cell_group == "Tumor cells")
## add label of SETD2 mutated or wt
### add readable aliquot id
plot_data_df$Id_Aliquot_WU <- mapvalues(x = plot_data_df$orig.ident, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
### map SETD2 mutation status
plot_data_df$Mut.SETD2 <- mapvalues(x = plot_data_df$Id_Aliquot_WU, from = bulkprofile_df$Aliquot_snRNA_WU, to = as.vector(bulkprofile_df$Mut.SETD2))
## filter
plot_data_df <- plot_data_df %>%
  filter(!is.na(Mut.SETD2)) %>%
  filter(Mut.SETD2 != Id_Aliquot_WU) %>%
  mutate(Is.Mut.SETD2 = ifelse(Mut.SETD2 == "None", "WT", "MT"))
## plot
p <- ggplot()
p <- p + geom_point(data = plot_data_df, 
                    mapping = aes(x = UMAP_1, y = UMAP_2, color = Is.Mut.SETD2),
                    alpha = 1, size = 0.05, shape = 16)
p <- p + ggtitle(label = "Tumor Cells By SETD2 Mutation Status") 
p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p <- p + theme(legend.position="bottom")
p
## save as png
file2write <- paste0(dir_out, "SETD2", ".png")
png(filename = file2write, width = 1000, height = 1100, res = 150)
print(p)
dev.off()

