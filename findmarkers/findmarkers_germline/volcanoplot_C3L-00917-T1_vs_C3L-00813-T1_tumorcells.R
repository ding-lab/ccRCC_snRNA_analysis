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
markers_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_germline/findallmarker_wilcox_germline_vhl_vs_2nd_somatic_nonsense_bycelltype_on_katmai/20200604.v1/VHL_Germline_vs_VHL_Somatic.FindMarkers.Wilcox.20200604.v1.tsv")
## input ppi table
highlight_degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_germline/examine_degs/overlap_C3L-00917-T1_vs_4_somatic_nonsense_VHL_mutants_tumorcells/20210303.v2/Consistent_DEGs_Wide.20210303.v2.tsv")

# preprocess------------------------------------------------------------------
## filter vhl interactome
genes_highlight_df <- highlight_degs_df %>%
  filter(is_gene_in_ora_pathways & (over3rdQu | below1stQu))
genes_highlight <- genes_highlight_df$deg_gene_symbol
## 
genegroup_highlight <- "P.adjusted<0.05, top fold change & enriched pathway member"

# make plot data  ---------------------------------------------------------
celltype_plot <- "Tumor cells"
plot_data_df <- markers_df %>%
  filter(Cell_type.shorter == celltype_plot) %>%
  # filter(p_val_adj < 0.05) %>%
  mutate(logPvalue = -log10(p_val_adj))
max_logp <- max(plot_data_df$logPvalue[!is.infinite(plot_data_df$logPvalue)])
## cat value
plot_data_df <- plot_data_df %>%
  mutate(y_plot = ifelse(is.infinite(logPvalue), max_logp, logPvalue)) %>%
  mutate(genegroup = ifelse(p_val_adj > 0.05, "P.adjusted>0.05",
                            ifelse(deg_gene_symbol %in% genes_highlight, genegroup_highlight, "P.adjusted<0.05")))
table(plot_data_df$genegroup)
# make plot dependencies --------------------------------------------------
colors_genegroup <- c("grey50", "black", brewer.pal(n = 3, name = "Set1")[2])
names(colors_genegroup) <- c("P.adjusted>0.05", "P.adjusted<0.05", genegroup_highlight)

# plot --------------------------------------------------------------------
## make color for highlighted genes
p <- ggplot()
p <- p + geom_vline(xintercept = 0, linetype = 2)
p <- p + geom_point(data = subset(plot_data_df, genegroup == "P.adjusted>0.05"), mapping = aes(x = avg_logFC, y = y_plot, color = genegroup), alpha = 1, shape = 16, size = 0.7)
p <- p + geom_point(data = subset(plot_data_df, genegroup == "P.adjusted<0.05"), mapping = aes(x = avg_logFC, y = y_plot, color = genegroup), alpha = 0.6, shape = 16, size = 0.7)
p <- p + geom_point(data = subset(plot_data_df, genegroup == genegroup_highlight), mapping = aes(x = avg_logFC, y = y_plot, color = genegroup), alpha = 1, shape = 16, size = 1.5)
p <- p + scale_color_manual(values = colors_genegroup)
p <- p + geom_text_repel(data = subset(plot_data_df, genegroup == genegroup_highlight),
                         mapping = aes(x = avg_logFC, y = y_plot, label = deg_gene_symbol), color = colors_genegroup[genegroup_highlight], 
                         force = 2, box.padding = unit(0.5, "lines"), fontface = "bold.italic", fontsize = 20)
p <- p + xlim(c(-1,1))
p <- p + ylim(c(0, max_logp*1.1))
p <- p + ggtitle(label = paste0("Differentially Expressed Genes in ", celltype_plot), 
                 subtitle = paste0("Germline-VHL-Mutated ", "(", "CPT0023690004", ")",
                                   " vs\nSomatic-VHL-Mutated ", "(", "CPT0015810004", ")"))
p <- p + xlab("ln(Fold Change of Average Expression)")
p <- p + ylab("-log10(Adjusted P-value)")
p <- p + theme_classic(base_size = 20)
p <- p + theme(legend.position = "bottom")
p <- p + theme(title = element_text(size = 15))
p <- p + guides(colour = guide_legend(override.aes = list(size=3), nrow = 3))
p

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "volcanoplot.2rd_somatic_samples", ".pdf")
pdf(file = file2write, width = 6, height = 8)
print(p)
dev.off()

file2write <- paste0(dir_out, "volcanoplot.2rd_somatic_samples", ".png")
png(filename = file2write, width = 1000, height = 1200, res = 150)
print(p)
dev.off()
