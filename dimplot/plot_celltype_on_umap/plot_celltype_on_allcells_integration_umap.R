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

# plot for cell group----------------------------------------------------------
plot_data_df <- merge(integrated_umap_df,
                      barcode2celltype_df %>%
                        select(integrated_barcode, Cell_group),
                      by.x = c("barcode"), by.y = c("integrated_barcode"), all.x = T
)

p <- ggplot()
p <- p + geom_point(data = plot_data_df, 
                    mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_group),
                    alpha = 1, size = 0.05, shape = 16)
p <- p + scale_color_manual(values = cellgroup_colors)
p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p <- p + theme(legend.position="bottom")
p
## save as pdf
file2write <- paste0(dir_out, "cellgroup_on_umap.", run_id, ".pdf")
pdf(file = file2write, width = 8, height = 9, useDingbats = F)
print(p)
dev.off()
## save as png
file2write <- paste0(dir_out, "cellgroup_on_umap.", run_id, ".png")
png(filename = file2write, width = 1000, height = 1100, res = 150)
print(p)
dev.off()

# plot for Cell_type.shorter all cells----------------------------------------------------------
plot_data_df <- merge(integrated_umap_df,
                      barcode2celltype_df %>%
                        select(integrated_barcode, Cell_type.shorter),
                      by.x = c("barcode"), by.y = c("integrated_barcode"), all.x = T
)

p <- ggplot()
p <- p + geom_point(data = plot_data_df, 
                    mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_type.shorter),
                    alpha = 1, size = 0.05)
p <- p + scale_color_manual(values = celltype_shorter_colors)
p <- p + guides(colour = guide_legend(override.aes = list(size=3)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p
file2write <- paste0(dir_out, "celltype_shorter_on_umap.", run_id, ".png")
png(filename = file2write, width = 1600, height = 1000, res = 150)
print(p)
dev.off()

# plot for Cell_type.shorter immune cells----------------------------------------------------------
plot_data_df <- merge(integrated_umap_df,
                      barcode2celltype_df %>%
                        select(integrated_barcode, Cell_type.shorter, Cell_group),
                      by.x = c("barcode"), by.y = c("integrated_barcode"), all.x = T
)
plot_data_df <- plot_data_df %>%
  filter(Cell_group == "Immune")

p <- ggplot()
p <- p + geom_point(data = plot_data_df, 
                    mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_type.shorter),
                    alpha = 1, size = 0.05)
p <- p + scale_color_manual(values = celltype_shorter_colors)
p <- p + guides(colour = guide_legend(override.aes = list(size=3)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p
file2write <- paste0(dir_out, "immune.", "celltype1_on_umap.", run_id, ".png")
png(filename = file2write, width = 1600, height = 1000, res = 150)
print(p)
dev.off()

# plot for Cell_type.shorter myeloid cells----------------------------------------------------------
plot_data_df <- merge(integrated_umap_df,
                      barcode2celltype_df %>%
                        select(integrated_barcode, Cell_type.shorter, Cell_group, Most_Enriched_Cell_Type1),
                      by.x = c("barcode"), by.y = c("integrated_barcode"), all.x = T
)
plot_data_df <- plot_data_df %>%
  filter(Most_Enriched_Cell_Type1 == "Myleoid lineage immune cells")

p <- ggplot()
p <- p + geom_point(data = plot_data_df, 
                    mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_type.shorter),
                    alpha = 1, size = 0.05)
p <- p + scale_color_manual(values = immune_myeloid_colors)
p <- p + guides(colour = guide_legend(override.aes = list(size=3)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p
file2write <- paste0(dir_out, "myleoid.", "celltype_shorter_on_umap.", run_id, ".png")
png(filename = file2write, width = 1200, height = 1000, res = 150)
print(p)
dev.off()


# plot for Cell_type.shorter lymloid cells----------------------------------------------------------
plot_data_df <- merge(integrated_umap_df,
                      barcode2celltype_df %>%
                        select(integrated_barcode, Cell_type.shorter, Cell_group, Most_Enriched_Cell_Type1),
                      by.x = c("barcode"), by.y = c("integrated_barcode"), all.x = T
)
plot_data_df <- plot_data_df %>%
  filter(Most_Enriched_Cell_Type1 == "Lymphoid lineage immune cells")

p <- ggplot()
p <- p + geom_point(data = plot_data_df, 
                    mapping = aes(x = UMAP_1, y = UMAP_2, color = Cell_type.shorter),
                    alpha = 1, size = 0.05)
p <- p + scale_color_manual(values = immune_lymphoid_colors)
p <- p + guides(colour = guide_legend(override.aes = list(size=3)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p
file2write <- paste0(dir_out, "lymphoid.", "celltype_shorter_on_umap.", run_id, ".png")
png(filename = file2write, width = 1200, height = 1000, res = 150)
print(p)
dev.off()

# plot for Most_Enriched_Cell_Type1 immune cells----------------------------------------------------------
plot_data_df <- merge(integrated_umap_df,
                      barcode2celltype_df %>%
                        select(integrated_barcode, Cell_type.shorter, Cell_group, Most_Enriched_Cell_Type1),
                      by.x = c("barcode"), by.y = c("integrated_barcode"), all.x = T
)
plot_data_df <- plot_data_df %>%
  filter(Cell_group == "Immune")

p <- ggplot()
p <- p + geom_point(data = plot_data_df, 
                    mapping = aes(x = UMAP_1, y = UMAP_2, color = Most_Enriched_Cell_Type1),
                    alpha = 1, size = 0.05)
p <- p + scale_color_manual(values = immunecelltype1_colors)
p <- p + guides(colour = guide_legend(override.aes = list(size=3)))
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.text.x=element_blank(),
               axis.ticks.x=element_blank())
p <- p + theme(axis.text.y=element_blank(),
               axis.ticks.y=element_blank())
p
file2write <- paste0(dir_out, "immune.", "celltype1_on_umap.", run_id, ".png")
png(filename = file2write, width = 1400, height = 1000, res = 150)
print(p)
dev.off()


