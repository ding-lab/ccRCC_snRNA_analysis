# Yige Wu @WashU Sep 2020

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

# input dependencies ------------------------------------------------------
## input EMT related genes
emtgenes_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/write_emt_genes/20200915.v1/EMT_Genes.20200915.v1.tsv")
## input cell type markers
gene2celltype_df <- fread(data.table = F, input = "./Resources/Knowledge/Kidney_Markers/Gene2CellType_Tab.20200911.v1.tsv")
# gene2celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/combine_pt_with_emt_markers_all/20200920.v1/Kidney_Specific_EMT_Genes.20200920.v1.tsv")

# filter deg --------------------------------------------------------------
# ## input degs
# deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/findmarkers_C3N-01200_vim_high_vs_low/20200910.v1/findmarkers_wilcox_vimhigh_vs_low.logfcthreshold0.693147180559945.minpct0.1.mindiffpct0.1.tsv")
# deg_df <- deg_df %>%
#   filter(clusterid_group1 == 6)
deg_df <- fread(data.table = F, input = "./Resources/snRNA_Processed_Data/Differentially_Expressed_Genes/C3N-01200.Tumor_Segments.Merged.DEGs.monocle.branch2.tsv")

# filter the emt genes ----------------------------------------------------
genes_filter <- unique(c(emtgenes_df$hgnc_symbol[emtgenes_df$gene_function != "Unknown"], gene2celltype_df$Gene[gene2celltype_df$Cell_Type1 %in% c("Proximal tubule")]) )
genes_filter

# make data for plotting --------------------------------------------------
deg_df$gene_function <- mapvalues(x = deg_df$Gene, from = emtgenes_df$hgnc_symbol, to = as.vector(emtgenes_df$gene_function))
deg_df$Cell_Type1 <- mapvalues(x = deg_df$Gene, from = gene2celltype_df$Gene, to = as.vector(gene2celltype_df$Cell_Type1))

plot_data_df <- deg_df %>%
  # dplyr::rename(genesymbol = row_name) %>%
  dplyr::rename(genesymbol = Gene) %>%
  mutate(Log10p_val_adj = -log10(x = p_val_adj))  %>%
  mutate(avg_log2FC = avg_logFC/log(2)) %>%
  mutate(x_plot = avg_log2FC) %>%
  mutate(text_gene = ifelse((genesymbol %in% genes_filter) & p_val_adj < 0.05 & ((gene_function == "Mesenchymal" & x_plot > 0) | (Cell_Type1 %in% c("Proximal tubule") & x_plot < 0)), genesymbol, NA))
## cap y axis
y_cap <- max(plot_data_df$Log10p_val_adj[!is.infinite(plot_data_df$Log10p_val_adj)])
plot_data_df <- plot_data_df %>%
  mutate(y_plot = ifelse(Log10p_val_adj >= y_cap, y_cap, Log10p_val_adj))
## set y bottom threshold
y_bottom <- -log10(0.05)
## set x limits to distinguish colors
x_pos <- log(2)
x_neg <- -log(2)

# plot --------------------------------------------------------------------
## plot
p <- ggplot()
# p <- p + geom_hline(yintercept = y_bottom, linetype = 2)
p <- p + geom_point(data = subset(plot_data_df, x_plot >= 0), mapping = aes(x = x_plot, y = y_plot), alpha = 0.5, size = 0.5, color = "red")
p <- p + geom_point(data = subset(plot_data_df, x_plot <= 0), mapping = aes(x = x_plot, y = y_plot), alpha = 0.5, size = 0.5, color = "blue")
# p <- p + geom_point(data = subset(plot_data_df, x_plot > x_neg & x_plot < x_pos), mapping = aes(x = x_plot, y = y_plot), alpha = 0.5, size = 0.5, color = "grey50")
# p <- p + scale_color_manual(values = c("FDR<0.05 (up)" = "red", "FDR<0.05 (down)" = "blue", "FDR<0.05" = "black", "FDR>=0.05" = "grey80"))
# p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_gene)),
#                          mapping = aes(x = x_plot, y = y_plot, label = text_gene), color = "black", force = 3, alpha = 0.8, fontface = "bold")
p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_gene) & x_plot > 0),
                         mapping = aes(x = x_plot, y = y_plot, label = text_gene), color = "black", force = 3, alpha = 0.8, fontface = "bold", xlim = c(1.5, NA))
p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_gene) & x_plot < 0),
                         mapping = aes(x = x_plot, y = y_plot, label = text_gene), color = "black", force = 2, alpha = 0.8, fontface = "bold", xlim = c(NA, -1))
p <- p + theme_bw()
# p <- p + ggtitle(label = paste0("C3N-01200 ", "VIM-high transitional cells vs tumor cells"))
p <- p + xlab("log2(Fold-Change) (transitional cells vs tumor cells)")
p <- p + ylab("-Log10(P-value-adjusted)")
p

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "volcano.", "png")
png(file2write, width = 800, height = 700, res = 150)
print(p)
dev.off()
file2write <- paste0(dir_out, "volcano.", "pdf")
pdf(file2write, width = 6, height = 5, useDingbats = F)
print(p)
dev.off()


