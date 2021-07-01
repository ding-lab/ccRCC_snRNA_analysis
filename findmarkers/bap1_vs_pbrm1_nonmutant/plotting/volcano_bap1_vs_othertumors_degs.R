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
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/bap1_vs_pbrm1_nonmutant/summarize_degs/unite_BAP1_snRNA_bulkRNA_protein_DEGs/20210614.v1/BAP1_snRNA_DEGs.CNVcorrected.20210614.v1.tsv")

# set plotting parameters -------------------------------------------------
color_purple <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[4]
color_red <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[1]
color_blue <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[2]

# make data for plotting --------------------------------------------------
plot_data_df <- deg_df %>%
  filter(!is.na(FDR.cnvcorrected)) %>%
  mutate(log10FDR = -log10(FDR.cnvcorrected)) %>%
  mutate(foldchange_type = ifelse(FDR.cnvcorrected < 0.05, ifelse(avg_lo2FC.alltumorcells > 0  & Num_down == 0 & Num_sig_up >= 5, "consistently higher in BAP1-mutants",
                                  ifelse(avg_lo2FC.alltumorcells < 0  & Num_up == 0 & Num_sig_down >= 5, "consistently lower in BAP1-mutants", "insignificant")), "insignificant")) %>%
  mutate(size_plot = abs(Num_sig_up - Num_sig_down)/24) %>%
  mutate(x_plot = ifelse(avg_lo2FC.alltumorcells < -2.5, -2.5, ifelse(avg_lo2FC.alltumorcells > 2.5, 2.5,  avg_lo2FC.alltumorcells))) %>%
  mutate(y_plot = ifelse(log10FDR > 350, 350, log10FDR)) %>%
  arrange(desc(foldchange_type))

# plot all markers size not scaled--------------------------------------------------------------------
## plot
p <- ggplot()
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey70")
p <- p + geom_point(data = subset(plot_data_df, foldchange_type == "insignificant"), 
                    mapping = aes(x = x_plot, y = y_plot, color = foldchange_type), alpha = 0.5, shape = 16)
p <- p + geom_point(data = subset(plot_data_df, foldchange_type %in% c("consistently higher in BAP1-mutants", "consistently lower in BAP1-mutants")), 
                    mapping = aes(x = x_plot, y = y_plot,  color = foldchange_type), alpha = 0.8, shape = 16)
p <- p + scale_color_manual(values = c("consistently higher in BAP1-mutants" = color_red, 
                                       "consistently lower in BAP1-mutants" = color_blue, 
                                       "insignificant" = "grey50"))
p <- p + geom_text_repel(data = subset(plot_data_df, foldchange_type == "consistently higher in BAP1-mutants" & x_plot >= 1),
                         mapping = aes(x = x_plot, y = y_plot, label = genesymbol_deg),
                         color = "black", alpha = 1, size = 4.5, #fontface = "bold",
                         segment.size = 0.2, segment.alpha = 0.5, min.segment.length = 0,
                         xlim = c(0, NA), max.overlaps = Inf)
p <- p + geom_text_repel(data = subset(plot_data_df, (foldchange_type == "consistently lower in BAP1-mutants") & (x_plot <= -1)),
                         mapping = aes(x = x_plot, y = y_plot, label = genesymbol_deg),
                         color = "black", alpha = 1, size = 4.5, #fontface = "bold",
                         segment.size = 0.2, segment.alpha = 0.5, min.segment.length = 0,
                         xlim = c(NA, 0), max.overlaps = Inf, force = 5)
p <- p + scale_size_area(max_size = 4)
p <- p + theme_classic()
p <- p + xlab("Log2(fold change) of tumor-cell sn gene expression for\nBAP1-mutant tumor (n=10) vs. other tumors(n=19)")
p <- p + ylab("-Log10FDR")
p <- p + ylim(c(0, 450)) + xlim(c(-2.5, 2.5))
p <- p + guides(color = guide_legend(title = "Fold change type", title.position = "top", title.theme = element_text(size = 14), 
                                     nrow = 4, override.aes = aes(size = 3), label.theme = element_text(size = 14)))
p <- p + theme(axis.text = element_text(size = 14),
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
