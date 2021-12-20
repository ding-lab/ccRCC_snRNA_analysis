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
x_cap <- Inf
version_tmp <- paste0("Cap", x_cap)
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the average expression
exp_wide_df <- fread(data.table = F, input = "./Resources/Analysis_Results/average_expression/avgexp_sct_data_bycelltypew_epithelial_bysample_katmai/20210907.v1/35_aliquot_merged.avgexp.SCT.data.Cell_group_w_epithelialcelltypes20210907.v1.tsv")
## input meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210809.v1/meta_data.20210809.v1.tsv")
## input genes
gene_plot_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/overlap_tumor_vs_pt_DEGs_w_tumor_vs_other_DEGs/20210702.v1/ccRCC_markers.Surface.20210702.v1.tsv")
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups_35aliquots/20210802.v1/35Aliquot.Barcode2CellType.20210802.v1.tsv", data.table = F)

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
## count cells
cellcount_df <- barcode2celltype_df %>%
  select(orig.ident, Cell_group_w_epithelialcelltypes) %>%
  table() %>%
  as.data.frame()
cellcount_df$cell_group.columnname <- mapvalues(x = cellcount_df$Cell_group_w_epithelialcelltypes, from = cellgroup_label_df$cell_type13, to = as.vector(cellgroup_label_df$cell_type13.columnname))
cellcount_df <- cellcount_df %>%
  mutate(id_sample_cell_group = paste0(orig.ident, "_", cell_group.columnname)) %>%
  mutate(keep = (Freq >= 10))
ids_sample_cell_group_keep <- cellcount_df$id_sample_cell_group[cellcount_df$keep]

# identify genes to plot -------------------------------------------------
genes_plot <- c("CA9", "RARRES2")
# for (gene_plot in unique(gene_plot_df$Gene)) {
# make plot data ----------------------------------------------------------
plotdata_wide_df <- exp_wide_df %>%
  filter(V1 %in% genes_plot)
plotdata_df <- melt(data = plotdata_wide_df)
summary(plotdata_df$value)
plotdata_df <- plotdata_df %>%
  mutate(id_sample_cell_group = gsub(x = variable, pattern = "SCT\\.", replacement = "")) %>%
  mutate(aliquot = str_split_fixed(string = id_sample_cell_group, pattern = "_", n = 2)[,1]) %>%
  mutate(cell_group.columnname = str_split_fixed(string = id_sample_cell_group, pattern = "_", n = 2)[,2]) %>%
  # mutate(x_plot = value)
  mutate(x_plot = ifelse(V1 == "CP" & value > 10, 10, value)) %>%
  filter(id_sample_cell_group %in% ids_sample_cell_group_keep)

plotdata_df$cell_group <- mapvalues(x = plotdata_df$cell_group.columnname, 
                                    from = cellgroup_label_df$cell_type13.columnname,
                                    to = as.vector(cellgroup_label_df$cell_type13))
plotdata_df$id_sample <- mapvalues(x = plotdata_df$aliquot, 
                                   from = idmetadata_df$Aliquot.snRNA,
                                   to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
ids_samples_sorted <- plotdata_df %>%
  filter(cell_group == "Tumor cells") %>%
  filter(V1 == "CA9") %>%
  arrange(desc(x_plot))
ids_samples_sorted <- unique(ids_samples_sorted$id_sample)
ids_samples_sorted <- c(ids_samples_sorted[!grepl(x = ids_samples_sorted, pattern = "\\-N")], "C3L-00079-N", "C3N-00242-N", "C3N-01200-N", "C3L-00088-N")
## add cell count
plotdata_df <- plotdata_df %>%
  filter(!(cell_group %in% c("Unknown", "Immune others", "Normal epithelial cells")))
# plotdata_df$y_plot <- plotdata_df$id_sample
plotdata_df$y_plot <- factor(x = plotdata_df$id_sample, levels = rev(ids_samples_sorted))
plotdata_df$V1 <- factor(x = plotdata_df$V1, levels = genes_plot)

# plot --------------------------------------------------------------------
p <- ggplot(data = plotdata_df, mapping = aes(x = y_plot, y = x_plot, fill = cell_group, color = as.character(cell_group == "Tumor cells")))
p <- p + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), alpha = 0.7, dotsize = 1.5)
p <- p + scale_fill_manual(values = colors_cellgroup[unique(plotdata_df$cell_group)])
p <- p + guides(fill = guide_legend(override.aes = list(size = 5)))
p <- p + scale_color_manual(values = c("TRUE" = "black", "FALSE" = NA))
p <- p + theme_classic(base_size = 12)
p <- p + coord_flip()
p <- p + facet_grid(cols = vars(V1), scales = "free_x")
p <- p + ylab("Normalized expression")
p <- p + theme(panel.grid.major.y = element_line(size=.1, color="grey50" ), strip.text = element_text(size = 15), strip.background = element_blank())
p <- p + theme(axis.text.y = element_text(size = 10), axis.title.y = element_blank())
p <- p + theme(axis.text.x = element_text(size = 10), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "last")))
# file2write <- paste0(dir_out, gene_plot, ".png")
# png(file2write, width = 800, height = 800, res = 150)
# print(p)
# dev.off()

file2write <- paste0(dir_out, paste0(genes_plot, collapse = "_"), ".pdf")
pdf(file2write, width = 11, height = 6.4, useDingbats = F)
print(p)
dev.off()
