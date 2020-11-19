# Yige Wu @WashU Feb 2020
## for plotting the fraction of immune/stroma/tumor cells/normal epithelial cells

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

# input dependencies ------------------------------------------------
## input barcodes to cell type
barcode2celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups/20201027.v1/31Aliquot.Barcode2CellType.20201027.v1.tsv")
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")

# make data for plotting --------------------------------------------------
## map cell groups
barcode2celltype_df <- barcode2celltype_df %>%
  mutate(Cell_group = Cell_group13)
## map easy ids
barcode2celltype_df$Aliquot_WU <- mapvalues(x = barcode2celltype_df$orig.ident, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
## sum the barcode fraction by cell group by sample
frac_barcodes_by_cellgroupaliquot_df <- barcode2celltype_df %>%
  group_by(Aliquot_WU, Cell_group) %>%
  summarise(Num_Barcodes_ByAliquotCellGroup = n())
num_barcodes <- barcode2celltype_df %>%
  group_by(Aliquot_WU) %>%
  summarise(Num_Barcodes_ByAliquot = n())
frac_barcodes_by_cellgroupaliquot_df$Num_Barcodes_ByAliquot <- mapvalues(x = frac_barcodes_by_cellgroupaliquot_df$Aliquot_WU, from = num_barcodes$Aliquot_WU, to = as.vector(num_barcodes$Num_Barcodes_ByAliquot))
frac_barcodes_by_cellgroupaliquot_df$Num_Barcodes_ByAliquot <- as.numeric(frac_barcodes_by_cellgroupaliquot_df$Num_Barcodes_ByAliquot)
frac_barcodes_by_cellgroupaliquot_df <- frac_barcodes_by_cellgroupaliquot_df %>%
  dplyr::mutate(Frac_Barcodes_ByAliquotCellGroup = Num_Barcodes_ByAliquotCellGroup/Num_Barcodes_ByAliquot)
## make data frame
plot_df <- frac_barcodes_by_cellgroupaliquot_df %>%
  mutate(Aliquot_Suffix = str_split_fixed(string = Aliquot_WU, pattern = "-", n = 3)[,3])
## label sample to tumor and normal
plot_df$Sample_Type <- mapvalues(x = plot_df$Aliquot_WU, from = idmetadata_df$Aliquot.snRNA.WU, to = idmetadata_df$Sample_Type)
plot_df$Case <- mapvalues(x = plot_df$Aliquot_WU, from = idmetadata_df$Aliquot.snRNA.WU, to = idmetadata_df$Case)
## map cell group
plot_df$Cell_group_facet <- mapvalues(x = plot_df$Cell_group, from = barcode2celltype_df$Cell_group, to = as.vector(barcode2celltype_df$Cell_group3))
## filter
plot_df <- plot_df %>%
  filter(Case %in% c("C3N-00733", "C3L-00416", "C3L-00088", "C3N-01200")) %>%
  filter(!(Cell_group %in% c("Unknown")))
## sort
plot_df$Case <- factor(x = plot_df$Case, levels = c("C3N-00733", "C3L-00416", "C3L-00088", "C3N-01200"))
plot_df$Cell_group_facet <- factor(x = plot_df$Cell_group_facet, levels = c("Nephron_Epithelium", "Immune", "Stroma"))
### sum the barcode fraction by cell group
frac_barcodes_by_cellgroup_df <- barcode2celltype_df %>%
  group_by(Cell_group) %>%
  summarise(Num_Barcodes_ByCellGroup = n()) %>%
  arrange(desc(Num_Barcodes_ByCellGroup))
plot_df$Cell_group <- factor(x = plot_df$Cell_group, levels = rev(frac_barcodes_by_cellgroup_df$Cell_group))


# plot data frame ---------------------------------------------------------
p <- ggplot()
p <- p + geom_point(data = plot_df, mapping = aes(x = Aliquot_WU, y = Cell_group, size = Frac_Barcodes_ByAliquotCellGroup, color = Aliquot_Suffix), shape = 16)
p <- p + facet_grid(Cell_group_facet~ Case, scales = "free", space = "free", drop = T)
p <- p + scale_color_manual(values = colors_tumor_segments)
p <- p + scale_size_area()
p <- p + scale_x_discrete(breaks = plot_df$Aliquot_WU, labels = plot_df$Aliquot_Suffix)
p <- p + xlab("T1-T3: Tumor Segments;\nN: Adjacent Normal Tissue;") + ylab("Fraction of Cell Group")
p <- p + theme(strip.text.x = element_text(angle = 90, face = "bold"),
               strip.background.x = element_blank(),
               strip.text.y = element_text(angle = 0),
               strip.background.y = element_blank())
p <- p + theme(panel.spacing = unit(0, "lines"),
               panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"),
               panel.border = element_rect(color = "black", fill = NA))
p <- p + theme(legend.position = "right")
p <- p + guides(fill = guide_legend(byrow = F, ncol = 1))
p <- p + theme(axis.text.x = element_text(size = 10, angle = 0, vjust = 0.5))
p <- p + theme(axis.title.y = element_blank())
p <- p + theme(axis.ticks = element_blank())
p
file2write <- paste0(dir_out, "Cell_Group_Composition.", run_id, ".pdf")
pdf(file2write, width = 7.5, height = 3.4, useDingbats = F)
print(p)
dev.off()

