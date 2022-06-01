# Yige Wu @WashU May 2022

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "ggplot2",
  "ggrastr",
  "ggrepel"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    print(paste0(pkg_name_tmp, "is being installed!"))
    BiocManager::install(pkgs = pkg_name_tmp, update = F)
    install.packages(pkg_name_tmp, dependencies = T)
  }
  print(paste0(pkg_name_tmp, " is installed!"))
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input -------------------------------------------------------------------
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/rna/knockout_cell_lines/DEG_analysis/run_edgeR_DEG_2caki1_cp_vs_2caki1_NT/20220517.v1/Caki1_CP_vs_Caki1_NT.DEGs.20220517.v1.tsv")
deg_down_filtered_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/rna/knockout_cell_lines/filter_DEGs/filter_2caki1_cp_vs_2caki_nt_degs/20220524.v1/CP_shRNA_down_degs.potential_pathogenic.20220524.v1.tsv")
# deg_down_filtered_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/rna/knockout_cell_lines/filter_DEGs/filter_2caki1_cp_vs_2caki_nt_degs/20220524.v2/CP_shRNA_down_degs.potential_pathogenic.20220524.v2.tsv")

# set plotting parameters -------------------------------------------------
color_red <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[1]
color_blue <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[2]
fdr_sig_cutoff <- 0.05
# genes_highlight <- c("CP", "VEGFC", "VEGFA")
# genes_highlight <- c(paste0("CXCL", c(1, 2, 3, 5, 6, 8)),
#                      paste0("CCL", c(2, 20, 28)),
#                      paste0("IL", c("1A", "1B", "15", "18", "7")))
# genes_highlight <- str_split(string = "CXCL6/IL6/CSF2/CXCL5/CXCL1/TNFRSF10A/PF4V1/IL18/OSMR/CXCL8/LTBR/TNFSF10/INHBA/CXCL2/TNFRSF1B/TNFSF14/CX3CL1/CCL5/PRLR/IL7R/TNFRSF11B/IL1A/TNFRSF10D/IL7/CCL2/CCL28/IL1B/CLCF1/IL20RB/TNFRSF14/CD27/PDGFRB/TGFB2/LIF/HGF/IL12RB1/IL1R1/CSF2RA/TNFRSF12A/CD70/GDF5/CCL20/VEGFC/MET/CD40/IL15/CSF1/EGFR/TNFSF9/CXCL3/FAS/TGFBR2/IL6ST/TNFRSF9/PDGFC/TNFRSF1A/IL4R/PLEKHO2/TNFRSF10B/EDA2R", pattern = "\\/")[[1]]
genes_highlight <- str_split(string = "EREG/CXCL6/IL6/CSF2/CXCL5/CLEC2B/HRNR/CXCL1/PI3/SERPINB7/COLEC10/CLEC4E/FLG/PF4V1/FLG2/IL18/S100A16/S100A6/CXCL8/SERPINE1/ADAMTS16/TNFSF10/INHBA/ANXA13/MMP7/BTC/SERPINA1/CXCL2/SERPINB9/S100A3/TNFSF14/PLAU/CX3CL1/CCL5/CTSZ/CTSS/GDF15/AREG/TGM2/ADAM28/FGF1/BMP5/IL1A/LGALS1/IL7/GREM1/TGFA/CCL2/CCL28/ADAMTS12/IL1B/ANGPTL4/P3H2/ADAM12/CLCF1/FGF2/MMP14/ANXA1/S100A2/BDNF/LOXL2/TGFB2/SLPI/LIF/SFRP4/HGF/MMP24/SERPINB8/ADAMTSL1/LGALS3/P4HA3/ANXA2/MMP19/GDF5/PAPPA2/EPGN/PLOD2/CCL20/VEGFC/FGFBP1/PLXNB3/S100A13/IL34/P4HTM/SDC4/TIMP2/SEMA3C/IL15/S100A11/CTSD/CTSB/CSF1/NRG1/CD109/MUC12/LGALS8/MMP2/ANXA4/TNFSF9/FSTL3/S100A4/P4HA1/HBEGF/PAPPA/SEMA4B/CXCL3/P4HA2/ANXA5/FSTL1/SDC1/ADAMTS3/EBI3/ELFN2/EGLN1/ANXA3/LOX/CTSO/FGF5/PDGFC/CTSC/ADAM9/ANXA6/ANXA7/SERPINH1/BMP1/ADAM10/ANXA11/S100A10/ADAM15/SERPINB6/CTSA/PLOD1/ADAM17", pattern = "\\/")[[1]]
genes_highlight <- genes_highlight[genes_highlight %in% deg_down_filtered_df$external_gene_name[deg_down_filtered_df$FDR < 0.05]]
# genes_highlight <- c(genes_highlight, "CP", "VEGFA", "VIM", "VCAM1", "ANGPTL4", "COL5A2", "COL8A1", "COL28A1")
genes_highlight <- c("CP", "COL4A1", "OSMR", "TGM2", "VEGFA")

# make data for plotting --------------------------------------------------
text_up <- paste0("Up (", length(which(plot_data_df$FDR < 0.05 & plot_data_df$logFC > 0)), ")")
text_down <- paste0("Down (", length(which(plot_data_df$FDR < 0.05 & plot_data_df$logFC < 0)), ")")

plot_data_df <- deg_df %>%
  mutate(log10FDR = -log10(FDR)) %>%
  mutate(foldchange_type = ifelse(FDR < fdr_sig_cutoff, ifelse(logFC > 0, text_up, text_down), "insignificant"))
summary(plot_data_df$log10FDR)
summary(plot_data_df$logFC)
table(plot_data_df$foldchange_type)
plot_data_df <- plot_data_df %>%
  # mutate(x_plot = ifelse(logFC < -3, -3, ifelse(logFC > 3, 3,  logFC))) %>%
  mutate(x_plot = logFC) %>%
  mutate(y_plot = ifelse(log10FDR > 350, 350, log10FDR)) %>%
  arrange(desc(foldchange_type))
summary(plot_data_df$logFC[plot_data_df$FDR < 0.05 & plot_data_df$logFC > 0])
x_pos_cutoff <- min(plot_data_df$logFC[plot_data_df$FDR < 0.05 & plot_data_df$logFC > 0])
x_neg_cutoff <- max(plot_data_df$logFC[plot_data_df$FDR < 0.05 & plot_data_df$logFC < 0])
## make colors
colors_deg <- c(color_red, color_blue, "grey50")
names(colors_deg) <- c(text_up, text_down, "insignificant")
# plot all markers size not scaled--------------------------------------------------------------------
## plot
p <- ggplot()
p <- p + geom_vline(xintercept = x_pos_cutoff, linetype = 2, color = "grey70")
p <- p + geom_vline(xintercept = x_neg_cutoff, linetype = 2, color = "grey70")
p <- p + geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey70")
p <- p + geom_point_rast(data = subset(plot_data_df, foldchange_type == "insignificant"), 
                         mapping = aes(x = x_plot, y = y_plot, color = foldchange_type), alpha = 0.5, shape = 16, size = 0.5)
p <- p + geom_point_rast(data = subset(plot_data_df, foldchange_type != "insignificant"), 
                         mapping = aes(x = x_plot, y = y_plot,  color = foldchange_type), alpha = 0.5, shape = 16, size = 0.5)
p <- p + scale_color_manual(values = colors_deg)
p <- p + geom_text_repel(data = subset(plot_data_df, logFC > 0 & (external_gene_name %in% genes_highlight)),
                         mapping = aes(x = x_plot, y = y_plot, label = external_gene_name),
                         color = "black", alpha = 1, size = 4.5, #fontface = "bold",
                         segment.size = 0.2, segment.alpha = 0.5, #min.segment.length = 0,
                         xlim = c(0, NA), max.overlaps = Inf)
p <- p + geom_text_repel(data = subset(plot_data_df, (logFC < 0) & (external_gene_name %in% genes_highlight)),
                         mapping = aes(x = x_plot, y = y_plot, label = external_gene_name),
                         color = "black", alpha = 1, size = 4.5, #fontface = "bold",
                         segment.size = 0.2, segment.alpha = 0.5, #min.segment.length = 0,
                         xlim = c(NA, 0), max.overlaps = Inf, force = 5)
p <- p + theme_classic()
# p <- p + xlab("Log2(fold change) of tumor-cell sn gene expression for\nBAP1-mutant tumor vs. other tumors(n=19)")
p <- p + xlab("Log2(fold change)")
p <- p + ylab("-Log10FDR")
# p <- p + xlim(c(-3, 3))
p <- p + guides(color = guide_legend(title = paste0("DEG (", length(which(deg_df$FDR < fdr_sig_cutoff)), ")"), 
                                     title.position = "top", title.theme = element_text(size = 14), 
                                     nrow = 4, override.aes = aes(size = 3), label.theme = element_text(size = 14)))
p <- p + theme(axis.text = element_text(size = 14, color = "black"),
               axis.title = element_text(size = 14),
               # legend.position = "right", 
               legend.position = c(0.80, 0.25),
               legend.box = "horizontal")
p
file2write <- paste0(dir_out, paste0(head(genes_highlight, 20), collapse = "_"), ".pdf")
pdf(file2write, width = 5, height = 4, useDingbats = F)
print(p)
dev.off()
