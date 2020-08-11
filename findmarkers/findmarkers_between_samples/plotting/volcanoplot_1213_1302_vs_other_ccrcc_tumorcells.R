# Yige Wu @WashU Aug 2020

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
markers_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_between_samples/findmarkers_1213_1302_vs_others_tumorcells/20200805.v1/findallmarkers_wilcox_1213_1302_vs_others..logfcthreshold0.1.minpct0.1.mindiffpct0.1.tsv")
## input interesting genes
genes_highlight_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/write_markers_by_intratumorheterogeneity_types/20200504.v1/markergenes_by_intratumorheterogeneity_types.20200504.v1.tsv")

# plot with selected gene highlighted ----------------------------------------------------
plot_data_df <- markers_df %>%
  mutate(logPvalue = -log10(p_val_adj))
max_logp <- max(plot_data_df$logPvalue[!is.infinite(plot_data_df$logPvalue)])
## cat value
plot_data_df <- plot_data_df %>%
  mutate(y_plot = ifelse(is.infinite(logPvalue), max_logp, logPvalue)) %>%
  mutate(genegroup = ifelse(p_val_adj > 0.05, "P.adjusted>0.05",
                            ifelse(gene %in% c(as.vector(genes_highlight_df$gene_symbol)), "P.adjusted<0.05, Important Process", "P.adjusted<0.05")))
table(plot_data_df$genegroup)
## make plot dependencies
colors_genegroup <- c("grey80", "black", "red")
names(colors_genegroup) <- c("P.adjusted>0.05", "P.adjusted<0.05", "P.adjusted<0.05, Important Process")
p <- ggplot()
p <- p + geom_point(data = plot_data_df, mapping = aes(x = avg_logFC, y = y_plot, color = genegroup), alpha = 0.7, shape = 16, size = 0.7)
p <- p + geom_text_repel(data = subset(plot_data_df, genegroup == "P.adjusted<0.05, Important Process"),
                         mapping = aes(x = avg_logFC, y = y_plot, label = gene), color = "red", force = 2)
p <- p + scale_color_manual(values = colors_genegroup)
# p <- p + xlim(c(-1,1))
p <- p + ylim(c(0, max_logp))
p <- p + ggtitle(label = paste0("C3N-01213 & C3L-130 vs the Rest"), subtitle = paste0("Differentially Expressed Genes in Tumor Cells Only"))
p <- p + xlab("log(Fold Change of Average Expression)\n(C3N-01213 & C3L-1302 vs Other ccRCC Tumor Cells)")
p <- p + ylab("-log10(Adjusted P-value)")
p <- p + theme_bw()
p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
p <- p + theme(legend.position = "right")
p
file2write <- paste0(dir_out, "volcanoplot.", "1213_1302_vs_other_ccrcc_tumorcells.", "selectedgenes.", ".pdf")
pdf(file = file2write, width = 10, height = 6, useDingbats = F)
print(p)
dev.off()

# plot with the largest fold change gene highlighted ----------------------------------------------------
plot_data_df <- markers_df %>%
  mutate(logPvalue = -log10(p_val_adj))
max_logp <- max(plot_data_df$logPvalue[!is.infinite(plot_data_df$logPvalue)])*1.1
## cat value
x_high <- quantile(x = plot_data_df$avg_logFC, prob = 0.99)
x_high
x_low <- quantile(x = plot_data_df$avg_logFC, prob = 0.01)
x_low
plot_data_df <- plot_data_df %>%
  mutate(y_plot = ifelse(is.infinite(logPvalue), max_logp, logPvalue)) %>%
  mutate(genegroup = ifelse(p_val_adj > 0.05, "P.adjusted>0.05", "P.adjusted<0.05")) %>%
  mutate(text_gene = ifelse((avg_logFC >= x_high | avg_logFC <= x_low) & p_val_adj < 0.05, gene, NA)) %>%
  mutate(text_color = ifelse(avg_logFC > 0, "red", "blue"))

table(plot_data_df$genegroup)
p <- ggplot()
p <- p + geom_point(data = plot_data_df[plot_data_df$p_val_adj > 0.05,], mapping = aes(x = avg_logFC, y = y_plot), color = "grey80", alpha = 1, shape = 16, size = 0.7)
p <- p + geom_point(data = plot_data_df[plot_data_df$p_val_adj <= 0.05,], mapping = aes(x = avg_logFC, y = y_plot), color = "black", alpha = 0.7, shape = 16, size = 0.7)
# p <- p + geom_point(data = plot_data_df, mapping = aes(x = avg_logFC, y = y_plot, color = genegroup), alpha = 0.7, shape = 16, size = 0.7)
p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_gene) & avg_logFC > 0),
                         mapping = aes(x = avg_logFC, y = y_plot, label = text_gene), color = "red", force = 2, alpha = 0.8)
p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_gene) & avg_logFC <= 0),
                         mapping = aes(x = avg_logFC, y = y_plot, label = text_gene), color = "blue", force = 2)
# p <- p + scale_color_manual(values = colors_genegroup)
p <- p +scale_color_gradient2(midpoint=0, low="blue", mid="white",
                              high="red", space ="Lab" )
p <- p + ylim(c(0, max_logp))
p <- p + ggtitle(label = paste0("C3N-01213 & C3L-130 vs the Rest"), subtitle = paste0("Differentially Expressed Genes in Tumor Cells Only"))
p <- p + xlab("log(Fold Change of Average Expression)\n(C3N-01213 & C3L-1302 vs Other ccRCC Tumor Cells)")
p <- p + ylab("-log10(Adjusted P-value)")
p <- p + theme_bw()
p <- p + guides(colour = guide_legend(override.aes = list(size=5)))
p <- p + theme(legend.position = "right")
p
file2write <- paste0(dir_out, "volcanoplot.", "1213_1302_vs_other_ccrcc_tumorcells.", "highfoldchange.", ".pdf")
pdf(file = file2write, width = 10, height = 6, useDingbats = F)
print(p)
dev.off()
file2write <- paste0(dir_out, "volcanoplot.", "1213_1302_vs_other_ccrcc_tumorcells.", "highfoldchange.", ".png")
png(filename = file2write, width = 1500, height = 1000, res = 150)
print(p)
dev.off()
