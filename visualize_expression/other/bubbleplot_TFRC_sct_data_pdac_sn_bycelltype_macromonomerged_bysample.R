# Yige Wu @WashU Apr 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA//"
setwd(dir_base)
source("./ccRCC_snRNA_analysis//load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/plotting.R")

## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the average expression
exp_long_df <- fread(data.table = F, input = "./Resources/snRNA_Processed_Data/Others/0.expression_meta_PDAC.txt")

# identify genes to plot -------------------------------------------------
genes_plot <- "TFRC"
## process cell type labels
cellgroup_label_df <- data.frame(cell_type13 = c("B-cells", "CD4+ T-cells", "CD8+ T-cells", "DC", "Endothelial cells", "Fibroblasts", "Immune others", "Macrophages", "NK cells", "Normal epithelial cells", "Tumor cells", "Unknown"))
cellgroup_label_df <- cellgroup_label_df %>%
  mutate(cell_type13.columnname = gsub(x = cell_type13, pattern = "\\-|\\+| ", replacement = "."))
x_cap <- 1.5

# make plot data ----------------------------------------------------------
plotdata_df <- exp_long_df %>%
  mutate(sample_suffix = str_split_fixed(string = label, pattern = "_", n = 3)[,2]) %>%
  mutate(id_sample = paste0(case, "_", sample_suffix)) %>%
  mutate(cell_group = ifelse(cell_type %in% c("Macrophage", "Monocyte"), "Macrophage", cell_type)) %>%
  group_by(id_sample, cell_group) %>%
  summarise(value = mean(TFRC)) %>%
  mutate(x_plot = ifelse(value > x_cap, x_cap, value))
summary(plotdata_df$value)
table(plotdata_df$cell_group)
ids_samples_sorted <- plotdata_df %>%
  filter(cell_group == "Macrophage") %>%
  arrange(desc(x_plot))
ids_samples_sorted <- ids_samples_sorted$id_sample
plotdata_df$y_plot <- factor(x = plotdata_df$id_sample, levels = rev(ids_samples_sorted))
plotdata_df <- rbind(plotdata_df[plotdata_df$cell_group != "Macrophage",], plotdata_df[plotdata_df$cell_group == "Macrophage",])
## make colors
colors_cellgroup <- Polychrome::palette36.colors(n = length(unique(plotdata_df$cell_group)))
names(colors_cellgroup) <- unique(plotdata_df$cell_group)
# colors_cellgroup["Macrophage"] <- "#FF7F00"
# plot --------------------------------------------------------------------
p <- ggplot(data = plotdata_df, mapping = aes(x = y_plot, y = x_plot, fill = cell_group, color = as.character(cell_group == "Macrophage")))
p <- p + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), alpha = 0.7)
p <- p + scale_fill_manual(values = colors_cellgroup)
p <- p + scale_color_manual(values = c("TRUE" = "black", "FALSE" = NA))
p <- p + theme_classic(base_size = 12)
p <- p + coord_flip()
p <- p + ylab("Normalized expression")
p <- p + theme(panel.grid.major.y = element_line(size=.1, color="black" ))
p <- p + theme(axis.text.y = element_text(size = 8), axis.title.y = element_blank())
p <- p + theme(axis.text.x = element_text(size = 12), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "last")))
file2write <- paste0(dir_out, "PDAC_TFRC", ".png")
png(file2write, width = 1200, height = 1000, res = 150)
print(p)
dev.off()
file2write <- paste0(dir_out, "PDAC_TFRC", ".pdf")
pdf(file2write, width = 7.5, height = 7, useDingbats = F)
print(p)
dev.off()
