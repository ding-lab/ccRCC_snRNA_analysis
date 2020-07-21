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
barcode2celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_to_all_cells_with_patch/20200720.v1/31AliquotIntegration.Barcode2CellType.TumorManualCluster.20200720.v1.tsv")
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")

# make data frame for plotting --------------------------------------------
plot_data_df <- barcode2celltype_df %>%
  filter(Cell_group == "Normal epithelial cells") %>%
  group_by(Cell_type.detailed) %>%
  summarise(value = n()) %>%
  mutate(prop = value/sum(plot_data_df$value)*100) %>%
  arrange(desc(prop)) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop) %>%
  mutate(text_label = paste0(round(x = prop, digits = 1), "%"))
plot_data_df$Cell_type.detailed <- factor(x = plot_data_df$Cell_type.detailed, levels = rev(plot_data_df$Cell_type.detailed))
sum(plot_data_df$Num_Barcodes_ByCellType)
unique(plot_data_df$Cell_type.detailed)
# make pie chart ----------------------------------------------------------
p <- ggplot(data = plot_data_df, aes(x="", y = prop, fill=Cell_type.detailed))
p <- p + geom_bar(stat="identity", width=1, color = "white")
p <- p + coord_polar("y", start=0)
p <- p + scale_fill_manual(values = normal_epithelial_colors)
# p <- p + scale_fill_brewer(palette="Set1")
p <- p + theme_void()
p <- p + geom_text(aes(y = ypos, label = text_label), color = "black", size=6)
p <- p + ggtitle(label = paste0("Cell Type Distribution of ", sum(plot_data_df$value), " Normal Epithelial Cells"))
p <- p + theme(legend.position = "bottom")
p
file2write <- paste0(dir_out, "pie_normal_epithelial_cells.pdf")
pdf(file2write, height = 6, width = 6)
print(p)
dev.off()

