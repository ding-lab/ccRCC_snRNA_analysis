## reference: https://www.cell.com/cell-stem-cell/pdf/S1934-5909(19)30166-3.pdf

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
daps_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/bap1/overlap_degs/overlap_bap1_vs_nonmutant_enhancer_promoter_peaks_28samples_with_degs/20211011.v1/BAP1_vs_NonMutant_DAP_DEG_Merged.20211011.v1.tsv")

# make data for plotting --------------------------------------------------
cutoff_log2FC_snATAC <- 0
plot_data_df <- daps_df %>%
  filter(!is.na(FDR.snATAC.cnvcorrected)) %>%
  mutate(log10FDR = -log10(FDR.snATAC.cnvcorrected)) %>%
  mutate(foldchange_type = ifelse(FDR.snATAC.cnvcorrected < 0.05, ifelse(avg_log2FC.snATAC > cutoff_log2FC_snATAC, "consistently higher in BAP1-mutants",
                                                                        ifelse(avg_log2FC.snATAC < (-cutoff_log2FC_snATAC), "consistently lower in BAP1-mutants", "insignificant")), 
                                  "insignificant")) %>%
  # mutate(size_plot = abs(Num_sig_up - Num_sig_down)/24) %>%
  mutate(x_plot = ifelse(avg_log2FC.snATAC < -2, -2, ifelse(avg_log2FC.snATAC > 2, 2,  avg_log2FC.snATAC))) %>%
  mutate(y_plot = ifelse(log10FDR > 350, 350, log10FDR)) %>%
  select(peak, foldchange_type, x_plot, y_plot) %>%
  unique()

plot_text_df <- daps_df %>%
  filter(!is.na(FDR.snATAC.cnvcorrected)) %>%
  filter(Gene %in% c("DLC1", "CES3", "PTPRJ", "RAPGEF5", "EPHA6",
                     "DDIT4", ## mTOR
                     "PEBP1", ## EGFR signaling
                     "DUSP1", ## TNFA signaling
                     "ETS1", "HSF4", "FAAP20", "AXL", ## bottom 5
                     "DMXL1", "FABP6", "SCN4A", "PTPRU", "CROCC" ## top 5
                     )) %>%
  group_by(Gene) %>%
  top_n(wt = abs(avg_log2FC.snATAC), n = 1) %>%
  select(peak, Gene) %>%
  unique()

plot_data_df$Gene <- mapvalues(x = plot_data_df$peak, from = plot_text_df$peak, to = as.vector(plot_text_df$Gene))
plot_data_df$Gene[plot_data_df$Gene == plot_data_df$peak] <- NA

# plot all markers size not scaled--------------------------------------------------------------------
color_red <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[1]
color_blue <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[2]

## plot
p <- ggplot()
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey70")
p <- p + geom_point_rast(data = subset(plot_data_df, foldchange_type == "insignificant"), 
                         mapping = aes(x = x_plot, y = y_plot, color = foldchange_type), alpha = 0.5, shape = 16)
p <- p + geom_point_rast(data = subset(plot_data_df, foldchange_type %in% c("consistently higher in BAP1-mutants", "consistently lower in BAP1-mutants")), 
                         mapping = aes(x = x_plot, y = y_plot,  color = foldchange_type), alpha = 0.8, shape = 16)
p <- p + scale_color_manual(values = c("consistently higher in BAP1-mutants" = color_red, 
                                       "consistently lower in BAP1-mutants" = color_blue, 
                                       "insignificant" = "grey50"))
p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(Gene) & foldchange_type == "consistently higher in BAP1-mutants" & x_plot >= cutoff_log2FC_snATAC),
                         mapping = aes(x = x_plot, y = y_plot, label = Gene),
                         color = "black", alpha = 1, size = 4.5, #fontface = "bold",
                         segment.size = 0.4, segment.alpha = 0.6, min.segment.length = 0, box.padding = 0.5,
                         xlim = c(0.8, NA), max.overlaps = Inf)
p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(Gene) & foldchange_type == "consistently lower in BAP1-mutants" & x_plot <= -(cutoff_log2FC_snATAC)),
                         mapping = aes(x = x_plot, y = y_plot, label = Gene),
                         color = "black", alpha = 1, size = 4.5, #fontface = "bold",
                         segment.size = 0.4, segment.alpha = 0.6, min.segment.length = 0, box.padding = 0.5,
                         xlim = c(NA, -1), max.overlaps = Inf, force = 5)
p <- p + scale_size_area(max_size = 4)
p <- p + theme_classic()
# p <- p + xlab("Log2(fold change) of tumor-cell sn gene expression for\nBAP1-mutant tumor vs. other tumors(n=19)")
p <- p + xlab("Log2(fold change) of tumor-cell peak accessibility")
p <- p + ylab("-Log10FDR")
# p <- p + ylim(c(0, 450)) + xlim(c(-2.5, 2.5))
p <- p + ylim(c(0, 360)) + xlim(c(-2.75, 2.5))
p <- p + guides(color = guide_legend(title = "Fold change type", title.position = "top", title.theme = element_text(size = 14), 
                                     nrow = 4, override.aes = aes(size = 3), label.theme = element_text(size = 14)))
p <- p + theme(axis.text = element_text(size = 14, color = "black"),
               axis.title = element_text(size = 14),
               legend.position = "bottom", legend.box = "horizontal")
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
# p <- p + geom_point(data = subset(plot_data_df, foldchange_type == "insignificant"), 
#                     mapping = aes(x = x_plot, y = y_plot, color = foldchange_type), alpha = 0.5, shape = 16)
# p <- p + geom_point(data = subset(plot_data_df, foldchange_type %in% c("consistently higher in BAP1-mutants", "consistently lower in BAP1-mutants")), 
#                     mapping = aes(x = x_plot, y = y_plot,  color = foldchange_type), alpha = 0.8, shape = 16)
# p <- p + scale_color_manual(values = c("consistently higher in BAP1-mutants" = color_red, 
#                                        "consistently lower in BAP1-mutants" = color_blue, 
#                                        "insignificant" = "grey50"))
# p <- p + geom_text_repel(data = subset(plot_data_df, foldchange_type == "consistently higher in BAP1-mutants" & x_plot >= 1),
#                          mapping = aes(x = x_plot, y = y_plot, label = Gene),
#                          color = "black", alpha = 1, size = 4.5, #fontface = "bold",
#                          segment.size = 0.2, segment.alpha = 0.5, min.segment.length = 0,
#                          xlim = c(0, NA), max.overlaps = Inf)
# p <- p + geom_text_repel(data = subset(plot_data_df, (foldchange_type == "consistently lower in BAP1-mutants") & (x_plot <= -1.25)),
#                          mapping = aes(x = x_plot, y = y_plot, label = Gene),
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
# 
# # write output ------------------------------------------------------------
# # file2write <- paste0(dir_out, "volcano.", "png")
# # png(file2write, width = 800, height = 800, res = 150)
# # print(p)
# # dev.off()
# file2write <- paste0(dir_out, "volcano.", "pdf")
# pdf(file2write, width = 5.5, height = 6.5, useDingbats = F)
# print(p)
# dev.off()
# 
