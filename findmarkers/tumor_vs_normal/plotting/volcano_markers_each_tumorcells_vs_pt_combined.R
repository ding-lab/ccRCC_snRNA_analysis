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
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/findallmarker_wilcox_each_atac_tumor_vs_pt_on_katmai/20201130.v2/findallmarkers_wilcox_each_snatac_tumor_vs_pt.20201130.v2.tsv")
## input degs that were commonly up-regulated
# deg_common_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/filter_deg/filter_each_tumor_vs_pt_degs_bysnatactumorgroup_shared/20201202.v1/Top_DEGs_Tumorcells_vs_PT_ByTumorGroup.20201202.v1.tsv")
deg_common_df <- fread(data.table = F, input = "./Resources/snRNA_Processed_Data/Differentially_Expressed_Genes/All_tumor-up_TF_DEG_withMotifs_inPromoterDACR.20201227.tsv")

# make plot data ----------------------------------------------------------
## summarize deg rsults
unique(deg_df$easyid_tumor)
unique(deg_df$genesymbol_deg) %>% length()
deg_common_filtered <- deg_common_df %>%
  filter(category_byshared_label == "all-tumor-shared")
# filter(direction_shared == "up" & category_byshared_label == "all-tumor-shared")
deg_filtered_df <- deg_df %>%
  # filter(genesymbol_deg %in% deg_common_filtered$genesymbol_deg) %>%
  filter(p_val_adj < 0.05)

## make plot data
plot_data_df <- deg_filtered_df %>%
  mutate(Log10p_val_adj = -log10(x = p_val_adj)) %>%
  group_by(genesymbol_deg) %>%
  summarise(avg_logFC_mean = mean(avg_logFC), Log10p_val_adj_mean = mean(Log10p_val_adj[!is.infinite(Log10p_val_adj)]), 
            number_upFC = sum(avg_logFC>0), number_downFC = sum(avg_logFC<0)) %>%
  mutate(FC_score = ifelse(avg_logFC_mean > 0, number_upFC, number_downFC)) %>%
  mutate(category_deg = ifelse(number_upFC >=6 & number_downFC == 0, "up", 
                               ifelse(number_upFC ==0 & number_downFC >=6, "down", "other"))) %>%
  arrange(FC_score)

# set plotting parameters --------------------------------------------------------------------
## subset genes to those to highlight
# genes_highlight <- c("EGLN3", "SLC38A1", "PFKP", "GATM", "KRBA1", "COL23A1", "TNIP1", "SLC15A4", "MXI1", "HILPDA", "EVA1C", "MSC-AS1",
#                      "PCK1", "PLXNB1", "SLC5A2", "ANPEP", "FGFR3", "NSMF", "BIN1", "SLC47A2")
genes_highlight1 <- plot_data_df %>%
  filter(genesymbol_deg %in% deg_common_df$V1) %>%
  arrange(desc(x_plot)) %>%
  head(n = 12)
genes_highlight <- c(genes_highlight1$genesymbol_deg,
                     "PCK1", "PLXNB1", "SLC5A2", "FGFR3", "NSMF", "SLC47A2")
plot_data_hl_df <- plot_data_df %>%
  filter(genesymbol_deg %in% genes_highlight)
color_up <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
color_down <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[2]
## set plotting parameters
x_cap_low <- -1.5
x_cap_high <- 1.5
x_lims <- c(x_cap_low, x_cap_high)
plot_data_df <- plot_data_df %>%
  mutate(x_plot = ifelse(avg_logFC_mean < x_cap_low, x_cap_low, ifelse(avg_logFC_mean > x_cap_high, x_cap_high, avg_logFC_mean))) %>%
  mutate(y_plot = Log10p_val_adj_mean)

# plot all markers--------------------------------------------------------------------
## plot
p <- ggplot()
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey50")
# p <- p + geom_point(data = plot_data_df, mapping = aes(x = x_plot, y = y_plot, alpha = FC_score), size = 0.5, color = "grey50")
# p <- p + geom_point(data = plot_data_df, mapping = aes(x = x_plot, y = y_plot), alpha = 0.8, size = 0.5, color = "red")
p <- p + geom_point(data = subset(plot_data_df, category_deg == "other"), mapping = aes(x = x_plot, y = y_plot), alpha = 0.1, size = 0.5, color = "grey50", shape = 16)
p <- p + geom_point(data = subset(plot_data_df, category_deg == "up"), mapping = aes(x = x_plot, y = y_plot), alpha = 0.6, size = 0.5, color = color_up, shape = 16)
p <- p + geom_point(data = subset(plot_data_df, category_deg == "down"), mapping = aes(x = x_plot, y = y_plot), alpha = 0.6, size = 0.5, color = color_down, shape = 16)
# p <- p + scale_color_manual(values = c("FDR<0.05 (up)" = "red", "FDR<0.05 (down)" = "blue", "FDR<0.05" = "black", "FDR>=0.05" = "grey80"))
p <- p + geom_text_repel(data = subset(plot_data_df, genesymbol_deg %in% genes_highlight & x_plot > 0),
                         mapping = aes(x = x_plot, y = y_plot, label = genesymbol_deg), 
                         force = 1, fontface = "italic", segment.alpha = 0.5, segment.size = 0.3,  xlim = c(0,NA))
p <- p + geom_text_repel(data = subset(plot_data_df, genesymbol_deg %in% genes_highlight & x_plot < 0),
                         mapping = aes(x = x_plot, y = y_plot, label = genesymbol_deg), 
                         force = 3, fontface = "italic", segment.alpha = 0.5, segment.size = 0.3, xlim = c(NA, 0))
p <- p + scale_color_manual(values = c("up" = color_up, "down" = color_down))
p <- p + xlim(x_lims)
# p <- p + ylim(c(0, y_top_area))
p <- p + xlab("ln(Fold Change)")
p <- p + ylab("-Log10(P-value-adjusted)")
p <- p + theme(legend.position = "bottom")
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank(), axis.line = element_line(colour = "black"))
p <- p + theme(axis.text = element_text(size = 12))
p
file2write <- paste0(dir_out, "tumorcells_vs_pt_DEGs", ".pdf")
pdf(file2write, width = 4, height = 3.5, useDingbats = F)
print(p)
dev.off()

