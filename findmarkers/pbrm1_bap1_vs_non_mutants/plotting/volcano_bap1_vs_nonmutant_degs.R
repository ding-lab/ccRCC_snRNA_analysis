# Yige Wu @WashU Sep 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(ggrastr)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/pbrm1_bap1_vs_non_mutants/summarize_degs/unite_BAP1_vs_NonMutant_snRNA_bulkRNA_protein_DEGs/20210913.v1/BAP1_snRNA_DEGs.CNVcorrected.20210913.v1.tsv")

# set plotting parameters -------------------------------------------------
color_purple <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[4]
color_red <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[1]
color_blue <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[2]
genes_highlight <- c("DLC1", #paste0("ARHGAP", c(24, 28, 32, 42)),
                     "PTPRJ", "CDH16", "CPEB3", "NR6A1", "ZBTB16",
                     "CES3", "PDK4", "SERPINA1", "SLC5A1", "TGFBR3",
                     "RAPGEF5", "MAPK9", "EPHA6", "EFNA5")
x_cutoff <- 2
y_cutoff <- 350
# make data for plotting --------------------------------------------------
plot_data_df <- deg_df %>%
  filter(!is.na(FDR.snRNA.cnvcorrected)) %>%
  mutate(log10FDR = -log10(FDR.snRNA.cnvcorrected)) %>%
  mutate(foldchange_type = ifelse(FDR.snRNA.cnvcorrected < 0.05, ifelse(avg_log2FC.snRNA > 0  & Num_down.snRNA == 0 & Num_sig_up.snRNA >= 5, "consistently higher in BAP1-mutants",
                                                                        ifelse(avg_log2FC.snRNA < 0  & Num_up.snRNA == 0 & Num_sig_down.snRNA >= 5, "consistently lower in BAP1-mutants", "insignificant")), "insignificant")) %>%
  # mutate(size_plot = abs(Num_sig_up - Num_sig_down)/24) %>%
  mutate(x_plot = ifelse(avg_log2FC.snRNA < (-x_cutoff), -x_cutoff, ifelse(avg_log2FC.snRNA > x_cutoff, x_cutoff,  avg_log2FC.snRNA))) %>%
  mutate(y_plot = ifelse(log10FDR > y_cutoff, y_cutoff, log10FDR)) %>%
  mutate(label_plot = ifelse(genesymbol_deg %in% genes_highlight, genesymbol_deg, NA)) %>%
  arrange(desc(foldchange_type))
table(plot_data_df$foldchange_type)
# plot all markers, highlight specified genes--------------------------------------------------------------------
## plot
p <- ggplot(data = plot_data_df, mapping = aes(x = x_plot, y = y_plot,  color = foldchange_type, label = label_plot))
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey70")
p <- p + geom_point_rast(alpha = 0.8, shape = 16)
p <- p + scale_color_manual(values = c("consistently higher in BAP1-mutants" = color_red, 
                                       "consistently lower in BAP1-mutants" = color_blue, 
                                       "insignificant" = "grey80"))
p <- p + geom_text_repel(data = subset(plot_data_df, x_plot < 0), 
                         color = "black", alpha = 1, size = 6, fontface = "italic",
                         segment.size = 0.2, segment.alpha = 0.5, min.segment.length = 0, box.padding = 0.75,
                         xlim = c(-x_cutoff, -0.3), force = 4,
                         max.overlaps = Inf)
p <- p + geom_text_repel(data = subset(plot_data_df, x_plot > 0), 
                         color = "black", alpha = 1, size = 6, fontface = "italic",
                         segment.size = 0.2, segment.alpha = 0.5, min.segment.length = 0, box.padding = 0.75,
                         xlim = c(0.2, x_cutoff), force = 4,
                         max.overlaps = Inf)
p <- p + scale_size_area(max_size = 4)
p <- p + theme_classic()
# p <- p + xlab("Log2(fold change) of tumor-cell sn gene expression for\nBAP1-mutant tumor vs. other tumors(n=19)")
p <- p + xlab("Log2(fold change) of tumor-cell gene expression")
p <- p + ylab("-Log10FDR")
# p <- p + ylim(c(0, 450)) + xlim(c(-2.5, 2.5))
p <- p + ylim(c(0, y_cutoff)) + xlim(c(-x_cutoff, x_cutoff))
p <- p + guides(color = guide_legend(title = "Fold change type", title.position = "top", title.theme = element_text(size = 14), 
                                     nrow = 4, override.aes = aes(size = 3), label.theme = element_text(size = 14)))
p <- p + theme(axis.text = element_text(size = 15, color = "black"),
               axis.title = element_text(size = 15),
               legend.position = "bottom", legend.box = "horizontal")
p
# file2write <- paste0(dir_out, "volcano.", "png")
# png(file2write, width = 800, height = 800, res = 150)
# print(p)
# dev.off()
file2write <- paste0(dir_out, "volcano.raster.", "pdf")
pdf(file2write, width = 5.5, height = 6.5, useDingbats = F)
print(p)
dev.off()

# # plot all markers size not scaled--------------------------------------------------------------------
# ## plot
# p <- ggplot()
# p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey70")
# p <- p + geom_point_rast(data = subset(plot_data_df, foldchange_type == "insignificant"), 
#                     mapping = aes(x = x_plot, y = y_plot, color = foldchange_type), alpha = 0.5, shape = 16)
# p <- p + geom_point_rast(data = subset(plot_data_df, foldchange_type %in% c("consistently higher in BAP1-mutants", "consistently lower in BAP1-mutants")), 
#                     mapping = aes(x = x_plot, y = y_plot,  color = foldchange_type), alpha = 0.8, shape = 16)
# p <- p + scale_color_manual(values = c("consistently higher in BAP1-mutants" = color_red, 
#                                        "consistently lower in BAP1-mutants" = color_blue, 
#                                        "insignificant" = "grey50"))
# p <- p + geom_text_repel(data = subset(plot_data_df, foldchange_type == "consistently higher in BAP1-mutants" & x_plot >= 1),
#                          mapping = aes(x = x_plot, y = y_plot, label = genesymbol_deg),
#                          color = "black", alpha = 1, size = 4.5, #fontface = "bold",
#                          segment.size = 0.2, segment.alpha = 0.5, min.segment.length = 0,
#                          xlim = c(0, NA), max.overlaps = Inf)
# p <- p + geom_text_repel(data = subset(plot_data_df, (foldchange_type == "consistently lower in BAP1-mutants") & (x_plot <= -1.25)),
#                          mapping = aes(x = x_plot, y = y_plot, label = genesymbol_deg),
#                          color = "black", alpha = 1, size = 4.5, #fontface = "bold",
#                          segment.size = 0.2, segment.alpha = 0.5, min.segment.length = 0,
#                          xlim = c(NA, 0), max.overlaps = Inf, force = 5)
# p <- p + scale_size_area(max_size = 4)
# p <- p + theme_classic()
# # p <- p + xlab("Log2(fold change) of tumor-cell sn gene expression for\nBAP1-mutant tumor vs. other tumors(n=19)")
# p <- p + xlab("Log2(fold change) of tumor-cell sn gene expression")
# p <- p + ylab("-Log10FDR")
# # p <- p + ylim(c(0, 450)) + xlim(c(-2.5, 2.5))
# p <- p + ylim(c(0, 450)) + xlim(c(-2.75, 2.5))
# p <- p + guides(color = guide_legend(title = "Fold change type", title.position = "top", title.theme = element_text(size = 14), 
#                                      nrow = 4, override.aes = aes(size = 3), label.theme = element_text(size = 14)))
# p <- p + theme(axis.text = element_text(size = 14, color = "black"),
#                axis.title = element_text(size = 14),
#                legend.position = "bottom", legend.box = "horizontal")
# # file2write <- paste0(dir_out, "volcano.", "png")
# # png(file2write, width = 800, height = 800, res = 150)
# # print(p)
# # dev.off()
# file2write <- paste0(dir_out, "volcano.raster.", "pdf")
# pdf(file2write, width = 5.5, height = 6.5, useDingbats = F)
# print(p)
# dev.off()
# 
# # plot all markers size not scaled--------------------------------------------------------------------
# ## plot
# p <- ggplot()
# p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey70")
# p <- p + geom_point(data = subset(plot_data_df, foldchange_type == "insignificant"), 
#                     mapping = aes(x = x_plot, y = y_plot, color = foldchange_type), alpha = 0.5, shape = 16)
# p <- p + geom_point(data = subset(plot_data_df, foldchange_type %in% c("consistently higher in BAP1-mutants", "consistently lower in BAP1-mutants")), 
#                     mapping = aes(x = x_plot, y = y_plot,  color = foldchange_type), alpha = 0.8, shape = 16)
# p <- p + scale_color_manual(values = c("consistently higher in BAP1-mutants" = color_red, 
#                                        "consistently lower in BAP1-mutants" = color_blue, 
#                                        "insignificant" = "grey50"))
# p <- p + geom_text_repel(data = subset(plot_data_df, foldchange_type == "consistently higher in BAP1-mutants" & x_plot >= 1),
#                          mapping = aes(x = x_plot, y = y_plot, label = genesymbol_deg),
#                          color = "black", alpha = 1, size = 4.5, #fontface = "bold",
#                          segment.size = 0.2, segment.alpha = 0.5, min.segment.length = 0,
#                          xlim = c(0, NA), max.overlaps = Inf)
# p <- p + geom_text_repel(data = subset(plot_data_df, (foldchange_type == "consistently lower in BAP1-mutants") & (x_plot <= -1.25)),
#                          mapping = aes(x = x_plot, y = y_plot, label = genesymbol_deg),
#                          color = "black", alpha = 1, size = 4.5, #fontface = "bold",
#                          segment.size = 0.2, segment.alpha = 0.5, min.segment.length = 0,
#                          xlim = c(NA, 0), max.overlaps = Inf, force = 5)
# p <- p + scale_size_area(max_size = 4)
# p <- p + theme_classic()
# # p <- p + xlab("Log2(fold change) of tumor-cell sn gene expression for\nBAP1-mutant tumor vs. other tumors(n=19)")
# p <- p + xlab("Log2(fold change) of tumor-cell sn gene expression")
# p <- p + ylab("-Log10FDR")
# p <- p + ylim(c(0, 450)) + xlim(c(-2.5, 2.5))
# p <- p + guides(color = guide_legend(title = "Fold change type", title.position = "top", title.theme = element_text(size = 14), 
#                                      nrow = 4, override.aes = aes(size = 3), label.theme = element_text(size = 14)))
# p <- p + theme(axis.text = element_text(size = 14, color = "black"),
#                axis.title = element_text(size = 14),
#                legend.position = "bottom", legend.box = "horizontal")
# # file2write <- paste0(dir_out, "volcano.", "png")
# # png(file2write, width = 800, height = 800, res = 150)
# # print(p)
# # dev.off()
# file2write <- paste0(dir_out, "volcano.", "pdf")
# pdf(file2write, width = 5.5, height = 6.5, useDingbats = F)
# print(p)
# dev.off()
