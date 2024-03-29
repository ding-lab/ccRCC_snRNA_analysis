# Yige Wu @WashU Apr 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA//"
setwd(dir_base)
source("./ccRCC_snRNA_analysis//load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/plotting.R")

## set run id
x_cap <- Inf
version_tmp <- paste0("Cap", x_cap)
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the average expression
exp_wide_df <- fread(data.table = F, input = "./Resources/Analysis_Results/average_expression/averexp_BRCA_scRNA_bycelltype_bysample_katmai/20210708.v1/BRCA.avgexp.SCT.data.cell_type.byPiece_ID.20210708.v1.tsv")

# preprocess --------------------------------------------------------------
## make colors
plotdata_wide_df <- exp_wide_df %>%
  filter(V1 %in% "CA9")
plotdata_df <- melt(data = plotdata_wide_df)
summary(plotdata_df$value)
plotdata_df <- plotdata_df %>%
  mutate(id_sample_cell_group = gsub(x = variable, pattern = "SCT\\.", replacement = "")) %>%
  mutate(case = str_split_fixed(string = id_sample_cell_group, pattern = "_", n = 3)[,1]) %>%
  mutate(sample_suffix = str_split_fixed(string = id_sample_cell_group, pattern = "_", n = 3)[,2]) %>%
  mutate(cell_group = str_split_fixed(string = id_sample_cell_group, pattern = "_", n = 3)[,3]) %>%
  # mutate(x_plot = value)
  mutate(x_plot = ifelse(value > x_cap, x_cap, value)) %>%
  mutate(id_sample = paste0(case, "_", sample_suffix))
colors_cellgroup <- Polychrome::palette36.colors(n = length(unique(plotdata_df$cell_group)))
names(colors_cellgroup) <- unique(plotdata_df$cell_group)

# identify genes to plot -------------------------------------------------
for (gene_plot in c("TF", "TFR2")) {
# for (gene_plot in c("TFRC", "HIF1A", "EPAS1", "SLC11A1", "SLC40A1")) {
  # make plot data ----------------------------------------------------------
  plotdata_wide_df <- exp_wide_df %>%
    filter(V1 %in% gene_plot)
  plotdata_df <- melt(data = plotdata_wide_df)
  summary(plotdata_df$value)
  plotdata_df <- plotdata_df %>%
    mutate(id_sample_cell_group = gsub(x = variable, pattern = "SCT\\.", replacement = "")) %>%
    mutate(case = str_split_fixed(string = id_sample_cell_group, pattern = "_", n = 3)[,1]) %>%
    mutate(sample_suffix = str_split_fixed(string = id_sample_cell_group, pattern = "_", n = 3)[,2]) %>%
    mutate(cell_group = str_split_fixed(string = id_sample_cell_group, pattern = "_", n = 3)[,3]) %>%
    # mutate(x_plot = value)
    mutate(x_plot = ifelse(value > x_cap, x_cap, value)) %>%
    mutate(id_sample = paste0(case, "_", sample_suffix))
  # plotdata_df$y_plot <- plotdata_df$id_sample
  
  ## sort by expression
  ids_samples_sorted_df <- plotdata_df %>%
    filter(cell_group == "Mono_Macro") %>%
    arrange(desc(x_plot))
  ids_samples_sorted <- unique(plotdata_df$id_sample)
  ids_samples_sorted <- c(ids_samples_sorted_df$id_sample, ids_samples_sorted[!(ids_samples_sorted %in% ids_samples_sorted_df$id_sample)])
  plotdata_df$y_plot <- factor(x = plotdata_df$id_sample, levels = rev(ids_samples_sorted))
  
  # plot --------------------------------------------------------------------
  p <- ggplot(data = plotdata_df, mapping = aes(x = y_plot, y = x_plot, fill = cell_group, color = as.character(cell_group == "Mono_Macro")))
  p <- p + geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(0.8), alpha = 0.7)
  p <- p + scale_fill_manual(values = colors_cellgroup[unique(plotdata_df$cell_group)])
  p <- p + scale_color_manual(values = c("TRUE" = "black", "FALSE" = NA))
  p <- p + theme_classic(base_size = 12)
  p <- p + coord_flip()
  p <- p + ylab("Normalized expression")
  p <- p + theme(panel.grid.major.y = element_line(size=.1, color="black" ))
  p <- p + theme(axis.text.y = element_text(size = 8), axis.title.y = element_blank())
  p <- p + theme(axis.text.x = element_text(size = 12), axis.line.x = element_line(arrow = grid::arrow(length = unit(0.3, "cm"), ends = "last")))
  p <- p + ggtitle(label = paste0(gene_plot, " sn Expression in BRCA"))
  # file2write <- paste0(dir_out, "BRCA_", gene_plot, ".png")
  # png(file2write, width = 1000, height = 1200, res = 150)
  # print(p)
  # dev.off()

  file2write <- paste0(dir_out, "BRCA_", gene_plot, ".pdf")
  pdf(file2write, width = 7, height = 5.5, useDingbats = F)
  print(p)
  dev.off()
}
