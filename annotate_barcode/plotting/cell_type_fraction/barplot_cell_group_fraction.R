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

# specify Cell group to plot ----------------------------------------------
var_cellgroup <- "Cell_group13"
colors_cellgroup <- colors_cellgroup13

# make data for plotting --------------------------------------------------
barcode2celltype_df$Aliquot_WU <- mapvalues(x = barcode2celltype_df$orig.ident, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
barcode2celltype_df[, "Cell_group"] <- barcode2celltype_df[, var_cellgroup]
## sum the barcode fraction by cell group
frac_barcodes_by_cellgroup <- barcode2celltype_df %>%
  group_by(Aliquot_WU, Cell_group) %>%
  summarise(Num_Barcodes_ByAliquotCellGroup = n())
num_barcodes <- barcode2celltype_df %>%
  group_by(Aliquot_WU) %>%
  summarise(Num_Barcodes_ByAliquot = n())
frac_barcodes_by_cellgroup$Num_Barcodes_ByAliquot <- mapvalues(x = frac_barcodes_by_cellgroup$Aliquot_WU, from = num_barcodes$Aliquot_WU, to = as.vector(num_barcodes$Num_Barcodes_ByAliquot))
frac_barcodes_by_cellgroup$Num_Barcodes_ByAliquot <- as.numeric(frac_barcodes_by_cellgroup$Num_Barcodes_ByAliquot)
frac_barcodes_by_cellgroup <- frac_barcodes_by_cellgroup %>%
  dplyr::mutate(Frac_Barcodes_ByAliquotCellGroup = Num_Barcodes_ByAliquotCellGroup/Num_Barcodes_ByAliquot)
## make data frame
plot_df <- frac_barcodes_by_cellgroup %>%
  mutate(Aliquot_Suffix = str_split_fixed(string = Aliquot_WU, pattern = "-", n = 3)[,3])
## label sample to tumor and normal
plot_df$Sample_Type <- mapvalues(x = plot_df$Aliquot_WU, from = idmetadata_df$Aliquot.snRNA.WU, to = idmetadata_df$Sample_Type)
plot_df$Case <- mapvalues(x = plot_df$Aliquot_WU, from = idmetadata_df$Aliquot.snRNA.WU, to = idmetadata_df$Case)

# sort cases by sample numbers and then normal epithelial cells--------------------------------------------------------------
case_sorted_df <- idmetadata_df %>%
  filter(snRNA_available) %>%
  select(Case) %>%
  table() %>%
  as.data.frame() %>%
  rename(Case = '.')
frac_normalepithelialcells_by_case <- plot_df %>%
  filter(Cell_group == "Normal epithelial cells") %>%
  group_by(Case) %>%
  summarise(Fraction_NormalEpiCell = sum(Frac_Barcodes_ByAliquotCellGroup))
case_sorted_df <- merge(case_sorted_df, frac_normalepithelialcells_by_case, by = c("Case"), all.x = T)
case_sorted_df$Fraction_NormalEpiCell[is.na(case_sorted_df$Fraction_NormalEpiCell)] <- 0
case_sorted_df <- case_sorted_df %>%
  arrange(Freq, Fraction_NormalEpiCell)
## sort
plot_df$Case <- factor(x = plot_df$Case, levels = case_sorted_df$Case)

# plot data frame ---------------------------------------------------------
p <- ggplot()
p <- p + geom_col(data = plot_df, mapping = aes(x = Aliquot_WU, y = Frac_Barcodes_ByAliquotCellGroup, fill = Cell_group), position = "stack")
p <- p + facet_grid(.~ Case, scales = "free", space = "free", drop = T)
p <- p + scale_fill_manual(values = colors_cellgroup)
p <- p + xlab("T1-T3: Tumor Segments;\nN: Adjacent Normal Tissue;") + ylab("Fraction of Cell Group")
p <- p + theme(strip.text.x = element_text(angle = 90, face = "bold"),
               strip.background = element_rect(color = "black", fill = "white"))
p <- p + scale_x_discrete(breaks = plot_df$Aliquot_WU, labels = plot_df$Aliquot_Suffix)
p <- p + theme(axis.text.x = element_text(size = 9, angle = 0, vjust = 0.5))
p <- p + theme(panel.spacing = unit(0, "lines"),
               panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"),
               panel.border = element_rect(color = "black", fill = NA, size = 1))
p <- p + theme(legend.position = "bottom")
p
file2write <- paste0(dir_out, "Cell_Group_Composition.", run_id, ".pdf")
pdf(file2write, width = 12, height = 5)
print(p)
dev.off()

