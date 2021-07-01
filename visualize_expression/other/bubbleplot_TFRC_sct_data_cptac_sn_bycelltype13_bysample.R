# Yige Wu @WashU Apr 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA//"
setwd(dir_base)
source("./ccRCC_snRNA_analysis//load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/plotting.R")
source("./ccRCC_snRNA_analysis/plotting_variables.R")

## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the average expression
exp_wide_df <- fread(data.table = F, input = "./Resources/Analysis_Results/average_expression/avgexp_sct_data_bycelltype13_bysample_katmai/20210701.v1/33_aliquot_merged.avgexp.SCT.data.bycelltype13_bysample.20210701.v1.tsv")
## input meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210423.v1/meta_data.20210423.v1.tsv")

# identify genes to plot -------------------------------------------------
genes_plot <- "TFRC"
## process cell type labels
cellgroup_label_df <- data.frame(cell_type13 = c("B-cells", "CD4+ T-cells", "CD8+ T-cells", "DC", "Endothelial cells", "Fibroblasts", "Immune others", "Macrophages", "NK cells", "Normal epithelial cells", "Tumor cells", "Unknown"))
cellgroup_label_df <- cellgroup_label_df %>%
  mutate(cell_type13.columnname = gsub(x = cell_type13, pattern = "\\-|\\+| ", replacement = "."))

# make plot data ----------------------------------------------------------
plotdata_wide_df <- exp_wide_df %>%
  filter(V1 %in% genes_plot)
plotdata_df <- melt(data = plotdata_wide_df)
summary(plotdata_df$value)
x_cap <- 5
plotdata_df <- plotdata_df %>%
  mutate(id_sample_cell_group = gsub(x = variable, pattern = "SCT\\.", replacement = "")) %>%
  mutate(aliquot = str_split_fixed(string = id_sample_cell_group, pattern = "_", n = 2)[,1]) %>%
  mutate(cell_group.columnname = str_split_fixed(string = id_sample_cell_group, pattern = "_", n = 2)[,2]) %>%
  mutate(x_plot = ifelse(value > x_cap, x_cap, value))

plotdata_df$cell_group <- mapvalues(x = plotdata_df$cell_group.columnname, 
                                    from = cellgroup_label_df$cell_type13.columnname,
                                    to = as.vector(cellgroup_label_df$cell_type13))
plotdata_df$id_sample <- mapvalues(x = plotdata_df$aliquot, 
                                    from = idmetadata_df$Aliquot.snRNA,
                                    to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
plotdata_df$y_plot <- plotdata_df$id_sample
plotdata_df <- plotdata_df %>%
  filter(!(cell_group %in% c("Unknown", "Immune others", "Normal epithelial cells")))
## make colors
colors_cellgroup <- c("#E7298A", "#1B9E77", "#7570B3","#000000", "#E69F00", "#56B4E9", "#F0E442", "#0072B2", "#D55E00","#CC79A7", "#B2DF8A", "#FB9A99","grey50")
names(colors_cellgroup) <- c("Tumor cells", "Normal epithelial cells", "Immune others", "B-cells", "CD4+ T-cells", "CD8+ T-cells", "Macrophages", "DC", "NK cells","Endothelial cells", "Myofibroblasts", "Fibroblasts","Unknown")

# plot --------------------------------------------------------------------
p <- ggplot(data = plotdata_df, mapping = aes(x = y_plot, y = x_plot, fill = cell_group, color = "white"))
p <- p + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), alpha = 0.7)
p <- p + scale_fill_manual(values = colors_cellgroup[unique(plotdata_df$cell_group)])
p <- p + scale_color_manual(values = c("white" = NA))
p <- p + theme_classic(base_size = 12)
p <- p + coord_flip()
p <- p + ylab("Normalized expression")
p <- p + theme(panel.grid.major.y = element_line(size=.1, color="black" ))
p <- p + scale_y_reverse()
p <- p + scale_x_discrete(position = "top")
p <- p + theme(legend.position = "left")
p <- p + theme(axis.text.y = element_text(size = 12), axis.title.y = element_blank())
p <- p + theme(axis.text.x = element_text(size = 12), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "first")))
p <- p + ggtitle(label = paste0("Human ccRCC sn Expression"))
file2write <- paste0(dir_out, "LigandGenes2", ".png")
png(file2write, width = 800, height = 800, res = 150)
print(p)
dev.off()


