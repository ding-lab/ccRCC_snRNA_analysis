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
version_tmp <- "Cap5"
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the average expression
exp_wide_df <- fread(data.table = F, input = "./Resources/Analysis_Results/average_expression/avgexp_sct_data_bycelltypew_epithelial_bysample_katmai/20210702.v1/33_aliquot_merged.avgexp.SCT.data.Cell_group_w_epithelialcelltypes20210702.v1.tsv")
## input meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210423.v1/meta_data.20210423.v1.tsv")
## input genes
gene_plot_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/overlap_tumor_vs_pt_DEGs_w_tumor_vs_other_DEGs/20210702.v1/ccRCC_markers.Surface.20210702.v1.tsv")
barcode2celltype_df <- fread(input = "./Data_Freezes/V2/snRNA/Cell_Type_Assignment/33Aliquot.Barcode2CellType.20210423.v1.tsv", data.table = F)

# preprocess --------------------------------------------------------------
## make colors
colors_tmp1 <- colors_cellgroup14[c("Tumor cells", "CD4+ T-cells", "CD8+ T-cells", "Macrophages", "NK cells", "DC", "Fibroblasts", "Myofibroblasts",  "B-cells")]
colors_tmp2 <- c(colors_cellgroup14[c("Normal epithelial cells", "EMT tumor cells", "Immune others")], "grey80",
                 Polychrome::palette36.colors(n = 36)[c("Vivid_Violet","Light_Olive_Brown", "Very_Light_Blue")], "grey80")
names(colors_tmp2) <- c("Proximal tubule", "Loop of Henle", "Distal convoluted tubule", 'Principle cells', 
                        "Intercalated cells", "Podocytes", "Endothelial cells", "Unknown")
colors_cellgroup <- c(colors_tmp1, colors_tmp2)
## process cell type labels
cellgroup_label_df <- data.frame(cell_type13 = c("B-cells", "CD4+ T-cells", "CD8+ T-cells", "DC", "Endothelial cells", "Fibroblasts", "Immune others", "Macrophages", "NK cells", 
                                                 "Normal epithelial cells", "Tumor cells", "Unknown",
                                                 "Proximal tubule", "Loop of Henle", "Distal convoluted tubule", 'Principle cells', "Intercalated cells", "Podocytes"))
cellgroup_label_df <- cellgroup_label_df %>%
  mutate(cell_type13.columnname = gsub(x = cell_type13, pattern = "\\-|\\+| ", replacement = "."))
x_cap <- 5

# identify genes to plot -------------------------------------------------
gene_plot <- gene_plot_df$Gene[1]
for (gene_plot in unique(gene_plot_df$Gene)) {
  # make plot data ----------------------------------------------------------
  plotdata_wide_df <- exp_wide_df %>%
    filter(V1 %in% gene_plot)
  plotdata_df <- melt(data = plotdata_wide_df)
  summary(plotdata_df$value)
  plotdata_df <- plotdata_df %>%
    mutate(id_sample_cell_group = gsub(x = variable, pattern = "SCT\\.", replacement = "")) %>%
    mutate(aliquot = str_split_fixed(string = id_sample_cell_group, pattern = "_", n = 2)[,1]) %>%
    mutate(cell_group.columnname = str_split_fixed(string = id_sample_cell_group, pattern = "_", n = 2)[,2]) %>%
    # mutate(x_plot = value)
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
  
  # plot --------------------------------------------------------------------
  p <- ggplot(data = plotdata_df, mapping = aes(x = y_plot, y = x_plot, fill = cell_group, color = "white"))
  p <- p + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), alpha = 0.7)
  p <- p + scale_fill_manual(values = colors_cellgroup[unique(plotdata_df$cell_group)])
  p <- p + scale_color_manual(values = c("white" = NA))
  p <- p + theme_classic(base_size = 12)
  p <- p + coord_flip()
  p <- p + ylab("Normalized expression")
  p <- p + theme(panel.grid.major.y = element_line(size=.1, color="black" ))
  p <- p + theme(axis.text.y = element_text(size = 12), axis.title.y = element_blank())
  p <- p + theme(axis.text.x = element_text(size = 12), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "last")))
  p <- p + ggtitle(label = paste0(gene_plot, " sn Expression"))
  file2write <- paste0(dir_out, gene_plot, ".png")
  png(file2write, width = 800, height = 800, res = 150)
  print(p)
  dev.off()
}
