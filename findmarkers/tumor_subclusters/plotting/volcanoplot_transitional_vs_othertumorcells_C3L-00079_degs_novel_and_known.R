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
## input degs
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/findmarkers_C3L-00079_transitional_vs_tumorcells_on_katmai/20201007.v1/findmarkers_wilcox_tumorlikecells_vs_tumorcells.logfcthreshold0.minpct0.1.mindiffpct0.1.tsv")
## input EMT related genes
gene2celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/combine_pt_with_emt_markers_all/20200920.v1/Kidney_Specific_EMT_Genes.20200920.v1.tsv")
emtgenes_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/write_emt_genes/20200915.v1/EMT_Genes.20200915.v1.tsv")
## input top pathway genes
ora_obj <- readRDS(file = "./Resources/Analysis_Results/pathway/clustreprofiler/clusterprofiler_deg_tumorlike_vs_tumorcells_C3L-00079/20200915.v1/ORA_Results.RDS")

# filter the genes to show ----------------------------------------------------
genes_mesenchymal <- emtgenes_df$hgnc_symbol[emtgenes_df$gene_function == "Mesenchymal"]
## 
ora_out_df <- as.data.frame(ora_obj)
genes_novel <- str_split(string = ora_obj@result$geneID[2], pattern = "\\/")[[1]]
genes_novel <- genes_novel[!(genes_novel %in% genes_mesenchymal)]
genes_filter <- c(genes_mesenchymal, genes_novel)
genes_filter

# set plotting parameters -------------------------------------------------
## set y bottom threshold
y_bottom <- -log10(0.05)
## set x limits to distinguish colors
x_pos <- 1
x_neg <- -1
## colors
color_red_deep <- RColorBrewer::brewer.pal(n = 6, name = "Paired")[6]
color_red_pale <- RColorBrewer::brewer.pal(n = 6, name = "Paired")[5]
color_blue_deep <- RColorBrewer::brewer.pal(n = 6, name = "Paired")[2]
color_blue_pale <- RColorBrewer::brewer.pal(n = 6, name = "Paired")[1]

# make data for plotting --------------------------------------------------
plot_data_df <- deg_df %>%
  dplyr::rename(genesymbol = row_name) %>%
  mutate(Log10p_val_adj = -log10(x = p_val_adj)) %>%
  mutate(avg_log2FC = avg_logFC/log(2)) %>%
  mutate(x_plot = avg_log2FC)
## cap y axis
y_cap <- max(plot_data_df$Log10p_val_adj[!is.infinite(plot_data_df$Log10p_val_adj)])
y_cap <- 200
plot_data_df <- plot_data_df %>%
  mutate(y_plot = ifelse(Log10p_val_adj >= y_cap, y_cap, Log10p_val_adj)) %>%
  mutate(text_gene = ifelse((y_plot >= y_bottom) & ((genesymbol %in% genes_mesenchymal) & (x_plot >= x_pos)) | ((genesymbol %in% genes_novel) & (x_plot <= x_neg | x_plot >= x_pos)), genesymbol, NA)) %>%
  mutate(font_face = ifelse(!is.na(text_gene) & text_gene %in% genes_mesenchymal, "bold.italic", "italic")) %>%
  mutate(font_color = ifelse(!is.na(text_gene) & text_gene %in% genes_mesenchymal, "red", "blue"))
  
# plot all markers--------------------------------------------------------------------
## plot
p <- ggplot()
p <- p + geom_vline(xintercept = x_pos, linetype = 2, alpha = 0.5)
p <- p + geom_vline(xintercept = x_neg, linetype = 2, alpha = 0.5)
p <- p + geom_point(data = subset(plot_data_df, y_plot < y_bottom), 
                    mapping = aes(x = x_plot, y = y_plot), alpha = 0.5, size = 1, color = "grey70")
p <- p + geom_point(data = subset(plot_data_df, x_plot < x_pos & x_plot > 0 & y_plot >= y_bottom), 
                    mapping = aes(x = x_plot, y = y_plot), alpha = 0.5, size = 1, color = color_red_pale)
p <- p + geom_point(data = subset(plot_data_df, x_plot >= x_pos & y_plot >= y_bottom),
                    mapping = aes(x = x_plot, y = y_plot), alpha = 0.5, size = 1, color = color_red_deep)
p <- p + geom_point(data = subset(plot_data_df, x_plot > x_neg & x_plot < 0 & y_plot >= y_bottom), 
                    mapping = aes(x = x_plot, y = y_plot), alpha = 0.5, size = 1, color = color_blue_pale)
p <- p + geom_point(data = subset(plot_data_df, x_plot <= x_neg & y_plot >= y_bottom), 
                    mapping = aes(x = x_plot, y = y_plot), alpha = 0.5, size = 1, color = color_blue_deep)
p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_gene) & x_plot > 0),
                         mapping = aes(x = x_plot, y = y_plot, label = text_gene, fontface = font_face), 
                         color = "black",
                         force = 5, segment.alpha = 0.5, size = 7, xlim = c(0.5, NA))
p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_gene) & x_plot < 0),
                         mapping = aes(x = x_plot, y = y_plot, label = text_gene, fontface = font_face), 
                         color = "black",
                         force = 5, segment.alpha = 0.5, size = 7, xlim = c(NA, -1))
p <- p + theme_bw()
# p <- p + xlim(c(-3, 3))
p <- p + xlab("log2(Fold-Change) (transitional cells vs tumor cells)")
p <- p + ylab("-Log10(P-value-adjusted)")
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank())
p <- p + theme(axis.title = element_text(size = 20))
p <- p + theme(legend.position = "none")
p <- p + theme(axis.text = element_text(size = 20))
p
file2write <- paste0(dir_out, "volcano.", "png")
png(file2write, width = 800, height = 600, res = 150)
print(p)
dev.off()
file2write <- paste0(dir_out, "volcano.", "pdf")
pdf(file2write, width = 9, height = 8, useDingbats = F)
print(p)
dev.off()


