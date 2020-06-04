# Yige Wu @WashU Jun 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(ggridges)
library(viridis)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input denpendencies -----------------------------------------------------
## input marker table
markers_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_germline/findallmarker_wilcox_germline_vhl_vs_1_somatic_nonsense_bycelltype_on_katmai/20200604.v1/VHL_Germline_vs_VHL_Somatic.FindMarkers.Wilcox.20200604.v1.tsv")
## input ppi table
ppi_pair_df <- fread(data.table = F, input = "./Resources/Databases/Protein_Protein_Interactions/protein_pair_table_v2.txt")

# get genes to highlight ------------------------------------------------------------------
## filter vhl interactome
genes_highlight_df <- ppi_pair_df %>%
  filter(GENE == "VHL")

# make plot data  ---------------------------------------------------------
celltype_plot <- "Tumor cells"
for (celltype_plot in c("Tumor cells")) {
# for (celltype_plot in c("Tumor cells", "Macrophages")) {
  max()
  plot_data_df <- markers_df %>%
    filter(Cell_type.shorter == celltype_plot) %>%
    # filter(p_val_adj < 0.05) %>%
    mutate(logPvalue = -log10(p_val_adj))
  max_logp <- max(plot_data_df$logPvalue[!is.infinite(plot_data_df$logPvalue)])
  ## cat value
  plot_data_df <- plot_data_df %>%
    mutate(y_plot = ifelse(is.infinite(logPvalue), max_logp, logPvalue)) %>%
    mutate(genegroup = ifelse(p_val_adj > 0.05, "P.adjusted>0.05",
                              ifelse(deg_gene_symbol %in% genes_highlight_df$SUB_GENE, "P.adjusted<0.05, Known to interact with VHL", "P.adjusted<0.05")))
  table(plot_data_df$genegroup)
  # make plot dependencies --------------------------------------------------
  colors_genegroup <- c("grey50", "black", "red")
  names(colors_genegroup) <- c("P.adjusted>0.05", "P.adjusted<0.05", "P.adjusted<0.05, Known to interact with VHL")
  
  # plot --------------------------------------------------------------------
  p <- ggplot()
  p <- p + geom_point(data = plot_data_df, mapping = aes(x = avg_logFC, y = y_plot, color = genegroup), alpha = 0.7, shape = 16, size = 0.7)
  p <- p + geom_text_repel(data = subset(plot_data_df, genegroup == "P.adjusted<0.05, Known to interact with VHL" & logPvalue > 100),
                           mapping = aes(x = avg_logFC, y = y_plot, label = deg_gene_symbol), color = "red", force = 2)
  p <- p + scale_color_manual(values = colors_genegroup)
  p <- p + xlim(c(-1,1))
  p <- p + ylim(c(0, max_logp))
  p <- p + ggtitle(label = paste0("Differentially Expressed Genes in ", celltype_plot), subtitle = paste0("Comparing Germline-VHL-Mutated Sample vs Somatic-VHL-Mutated Sample"))
  p <- p + xlab("log(Fold Change of Average Expression)\n(Germline-VHL-Mutated Sample vs Somatic-VHL-Mutated Sample)")
  p <- p + ylab("-log10(Adjusted P-value)")
  p <- p + theme_bw()
  p <- p + theme(legend.position = "bottom")
  p
  
  # save output -------------------------------------------------------------
  file2write <- paste0(dir_out, "volcanoplot.", celltype_plot, ".", "VHL_Germline", "_vs_", "VHL_Somatic.", run_id, ".pdf")
  pdf(file = file2write, width = 6, height = 6)
  print(p)
  dev.off()
  
  file2write <- paste0(dir_out, "volcanoplot.", celltype_plot, ".", "VHL_Germline", "_vs_", "VHL_Somatic", run_id,  ".png")
  png(filename = file2write, width = 1000, height = 1000, res = 150)
  print(p)
  dev.off()
}

for (celltype_plot in c("Macrophages")) {
  # for (celltype_plot in c("Tumor cells", "Macrophages")) {
  max()
  plot_data_df <- markers_df %>%
    filter(Cell_type.shorter == celltype_plot) %>%
    # filter(p_val_adj < 0.05) %>%
    mutate(logPvalue = -log10(p_val_adj))
  max_logp <- max(plot_data_df$logPvalue[!is.infinite(plot_data_df$logPvalue)])
  ## cat value
  plot_data_df <- plot_data_df %>%
    mutate(y_plot = ifelse(is.infinite(logPvalue), max_logp, logPvalue)) %>%
    mutate(genegroup = ifelse(p_val_adj > 0.05, "P.adjusted>0.05",
                              ifelse(deg_gene_symbol %in% genes_highlight_df$SUB_GENE, "P.adjusted<0.05, Known to interact with VHL", "P.adjusted<0.05")))
  table(plot_data_df$genegroup)
  # make plot dependencies --------------------------------------------------
  colors_genegroup <- c("grey50", "black", "red")
  names(colors_genegroup) <- c("P.adjusted>0.05", "P.adjusted<0.05", "P.adjusted<0.05, Known to interact with VHL")
  
  # plot --------------------------------------------------------------------
  p <- ggplot()
  p <- p + geom_point(data = plot_data_df, mapping = aes(x = avg_logFC, y = y_plot, color = genegroup), alpha = 0.7, shape = 16, size = 0.7)
  p <- p + geom_text_repel(data = subset(plot_data_df, genegroup == "P.adjusted<0.05, Known to interact with VHL"),
                           mapping = aes(x = avg_logFC, y = y_plot, label = deg_gene_symbol), color = "red", force = 2)
  p <- p + scale_color_manual(values = colors_genegroup)
  p <- p + xlim(c(-1,1))
  p <- p + ylim(c(0, max_logp))
  p <- p + ggtitle(label = paste0("Differentially Expressed Genes in ", celltype_plot), subtitle = paste0("Comparing Germline-VHL-Mutated Sample vs Somatic-VHL-Mutated Sample"))
  p <- p + xlab("log(Fold Change of Average Expression)\n(Germline-VHL-Mutated Sample vs Somatic-VHL-Mutated Sample)")
  p <- p + ylab("-log10(Adjusted P-value)")
  p <- p + theme_bw()
  p <- p + theme(legend.position = "bottom")
  p
  
  # save output -------------------------------------------------------------
  file2write <- paste0(dir_out, "volcanoplot.", celltype_plot, ".", "VHL_Germline", "_vs_", "VHL_Somatic.", run_id, ".pdf")
  pdf(file = file2write, width = 6, height = 6)
  print(p)
  dev.off()
  
  file2write <- paste0(dir_out, "volcanoplot.", celltype_plot, ".", "VHL_Germline", "_vs_", "VHL_Somatic", run_id,  ".png")
  png(filename = file2write, width = 1000, height = 1000, res = 150)
  print(p)
  dev.off()
}
