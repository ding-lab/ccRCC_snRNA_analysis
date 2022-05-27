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
genes_highlight <- str_split(string = "CXCL6/IL6/SPP1/CXCL1/NNMT/THBS2/PRRX1/CDH6/MYL9/CXCL8/VCAM1/SERPINE1/INHBA/LUM/CD44/TGFBI/NT5E/ANPEP/AREG/FMOD/TGM2/TFPI2/THBS1/ITGB3/COL6A3/TNFRSF11B/CAPG/IGFBP3/LGALS1/SPARC/GREM1/COL16A1/SPOCK1/ADAM12/IL32/MXRA5/COL5A2/OXTR/MYLK/FGF2/MMP14/DKK1/CCN2/BDNF/PDGFRB/LOXL2/TPM1/LAMC2/FN1/LAMA2/SFRP4/MATN2/GLIPR1/SAT1/ECM1/TNFRSF12A/PLOD2/ITGB1/VEGFC/TPM2/SDC4/CCN1/DAB2/IL15/PTX3/RGS4/VIM/QSOX1/MMP2/MEST/ITGAV/CD59/COL4A1/VCAN/EFEMP2/FSTL3/LRP1/BASP1/TAGLN/COL7A1/IGFBP4/DPYSL3/FBN1/PMEPA1/FAS/PLAUR/CALU/FSTL1/ITGA5/SDC1/COL6A2/ITGB5/COL4A2/EMP3/LOX/DST/CAP2/ENO2/MATN3/CALD1/LAMA1/SERPINH1/BMP1/ACTA2/COPA/PVR/TPM4/LAMA3/GEM/PLOD1/PPIB", pattern = "\\/")[[1]]
genes_highlight <- genes_highlight[genes_highlight %in% deg_down_filtered_df$external_gene_name[deg_down_filtered_df$FDR < 0.05]]
# genes_highlight <- c(genes_highlight, "CP", "VEGFA", "VIM", "VCAM1", "ANGPTL4", "COL5A2", "COL8A1", "COL28A1")
# make data for plotting --------------------------------------------------
# text_up <- paste0("Up (", len)
plot_data_df <- deg_df %>%
  mutate(log10FDR = -log10(FDR)) %>%
  mutate(foldchange_type = ifelse(FDR < fdr_sig_cutoff, ifelse(logFC > 0, "up", "down"), "insignificant"))
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

# plot all markers size not scaled--------------------------------------------------------------------
## plot
p <- ggplot()
p <- p + geom_vline(xintercept = x_pos_cutoff, linetype = 2, color = "grey70")
p <- p + geom_vline(xintercept = x_neg_cutoff, linetype = 2, color = "grey70")
p <- p + geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey70")
p <- p + geom_point_rast(data = subset(plot_data_df, foldchange_type == "insignificant"), 
                         mapping = aes(x = x_plot, y = y_plot, color = foldchange_type), alpha = 0.5, shape = 16)
p <- p + geom_point_rast(data = subset(plot_data_df, foldchange_type != "insignificant"), 
                         mapping = aes(x = x_plot, y = y_plot,  color = foldchange_type), alpha = 0.8, shape = 16)
p <- p + scale_color_manual(values = c("up" = color_red, 
                                       "down" = color_blue, 
                                       "insignificant" = "grey50"))
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
               legend.position = "right", legend.box = "horizontal")
p
file2write <- paste0(dir_out, paste0(head(genes_highlight, 20), collapse = "_"), ".pdf")
pdf(file2write, width = 6, height = 4, useDingbats = F)
print(p)
dev.off()
