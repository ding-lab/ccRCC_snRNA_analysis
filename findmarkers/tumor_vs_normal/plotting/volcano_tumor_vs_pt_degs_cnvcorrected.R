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
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
# deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/summarize_deg/unite_tumor_vs_normal_snRNA_individual_and_CNVcorrected_DEGs/20210608.v1/Tumor_vs_PT_DEGs.with.CNVcorrection.20210608.v1.tsv")
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/summarize_deg/unite_tumor_vs_normal_snRNA_individual_and_CNVcorrected_DEGs/20210824.v1/Tumor_vs_PT_DEGs.with.CNVcorrection.20210824.v1.tsv")

# set plotting parameters -------------------------------------------------
color_purple <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[4]
color_red <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[1]
color_blue <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[2]

# make data for plotting --------------------------------------------------
plot_data_df <- deg_df %>%
  filter(!is.na(FDR.CNVcorrected)) %>%
  mutate(log10FDR = -log10(FDR.CNVcorrected)) %>%
  mutate(foldchange_type = ifelse(FDR.CNVcorrected < 0.05, ifelse(avg_log2FC.allTumorcellsvsPT > 0  & Num_down == 0 & Num_sig_up >= 15, "consistently higher in ccRCC cells",
                                  ifelse(avg_log2FC.allTumorcellsvsPT < 0  & Num_up == 0 & Num_sig_down >= 15, "consistently lower in ccRCC cells", "insignificant")), "insignificant")) %>%
  # mutate(size_plot = abs(Num_sig_up - Num_sig_down)/24) %>%
  mutate(x_plot = ifelse(avg_log2FC.allTumorcellsvsPT < -3, -3, ifelse(avg_log2FC.allTumorcellsvsPT > 3, 3,  avg_log2FC.allTumorcellsvsPT))) %>%
  # mutate(x_plot = ifelse(avg_log2FC.allTumorcellsvsPT < -2.5, -2.5, ifelse(avg_log2FC.allTumorcellsvsPT > 2.5, 2.5,  avg_log2FC.allTumorcellsvsPT))) %>%
  mutate(y_plot = ifelse(log10FDR > 350, 350, log10FDR)) %>%
  arrange(desc(foldchange_type))

# plot all markers size not scaled--------------------------------------------------------------------
x_high <- 2.5
x_low <- -2.5
## plot
p <- ggplot()
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey70")
# p <- p + geom_vline(xintercept = x_high, linetype = 2, color = "grey70")
# p <- p + geom_vline(xintercept = x_low, linetype = 2, color = "grey70")
p <- p + geom_point(data = subset(plot_data_df, foldchange_type == "insignificant"), 
                    mapping = aes(x = x_plot, y = y_plot, color = foldchange_type), alpha = 0.5, shape = 16)
p <- p + geom_point(data = subset(plot_data_df, foldchange_type %in% c("consistently higher in ccRCC cells", "consistently lower in ccRCC cells")), 
                    mapping = aes(x = x_plot, y = y_plot,  color = foldchange_type), alpha = 0.8, shape = 16)
p <- p + scale_color_manual(values = c("consistently higher in ccRCC cells" = color_red, 
                                       "consistently lower in ccRCC cells" = color_blue, 
                                       "insignificant" = "grey50"))
p <- p + geom_text_repel(data = subset(plot_data_df, foldchange_type == "consistently higher in ccRCC cells" & x_plot >= x_high),
                         mapping = aes(x = x_plot, y = y_plot, label = genesymbol_deg),
                         color = "black", alpha = 1, size = 4.5, #fontface = "bold",
                         segment.size = 0.2, segment.alpha = 0.5, min.segment.length = 0,
                         xlim = c(0, NA), max.overlaps = Inf)
p <- p + geom_text_repel(data = subset(plot_data_df, (foldchange_type == "consistently lower in ccRCC cells") & (x_plot <= x_low)),
                         mapping = aes(x = x_plot, y = y_plot, label = genesymbol_deg),
                         color = "black", alpha = 1, size = 4.5, #fontface = "bold",
                         segment.size = 0.2, segment.alpha = 0.5, min.segment.length = 0,
                         xlim = c(NA, 0), max.overlaps = Inf, force = 5)
p <- p + scale_size_area(max_size = 4)
p <- p + theme_classic()
p <- p + xlab("Log2(fold change) of sn gene expression for\nccRCC cells vs. PT cells")
p <- p + ylab("-Log10FDR")
p <- p + ylim(c(0, 450)) + xlim(c(-3, 3))
p <- p + guides(color = guide_legend(title = "Fold change type", title.position = "top", title.theme = element_text(size = 14), 
                                     nrow = 4, override.aes = aes(size = 3), label.theme = element_text(size = 14)))
p <- p + theme(axis.text = element_text(size = 14, color = "black"),
               axis.title = element_text(size = 14),
               legend.position = "bottom", legend.box = "horizontal")
# file2write <- paste0(dir_out, "volcano.", "png")
# png(file2write, width = 800, height = 800, res = 150)
# print(p)
# dev.off()
file2write <- paste0(dir_out, "volcano.", "pdf")
pdf(file2write, width = 5.5, height = 6.5, useDingbats = F)
print(p)
dev.off()

p <- p + theme(legend.position = "none")
file2write <- paste0(dir_out, "volcano.nolgend.", "pdf")
pdf(file2write, width = 5.5, height = 5, useDingbats = F)
print(p)
dev.off()
