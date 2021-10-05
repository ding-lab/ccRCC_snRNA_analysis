# Yige Wu @WashU Sep 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(ggrastr)

## set run id
version_tmp <- 3
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input degs
dam_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Enriched_Motifs/EMT/Score_difference.EpithelialSelectedClusters_vs_Mesenchymal.20210924.tsv")
tf2deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/emt/filter_motifs/filter_motifs_for_EMT_vs_selectedEpithelial_peaks_for_degs/20210928.v1/MotifsTF2Foldchange.snRNA.20210928.v1.tsv")

# set plotting parameters -------------------------------------------------
## set y bottom threshold
y_bottom <- -log10(0.05)
## colors
color_right_deep <- RColorBrewer::brewer.pal(n = 12, name = "Set1")[1]
color_left_deep <- RColorBrewer::brewer.pal(n = 6, name = "Set1")[2]

# make data for plotting --------------------------------------------------
plot_data_df <- dam_df %>%
  rename(TF = TF_Name) %>%
  mutate(Is_FOSJUN = grepl(pattern = "FOS|JUN", x = TF) | (TF %in% c("ATF2", "ATF7", "ATF3", "ATF4", "BACH2", "BACH1", "BACH2(var.2)", "NFE2", "MAF::NFE2", "NFE2L1")) | (TF %in% c("JDP2", "JDP2(var.2)"))) %>%
  mutate(Log10p_val_adj = -log10(x = FDR)) %>%
  mutate(diff_2_vs_1 = (-diff)) %>%
  mutate(x_plot = ifelse((-diff) < -1.5, -1.5,
                         ifelse((-diff) > 1.5, 1.5, -diff))) %>%
  arrange(desc(diff_2_vs_1))
plot_data_df$avg_log2FC.tf.snRNA <- mapvalues(x = plot_data_df$TF, from = tf2deg_df$TF_name, to = as.vector(tf2deg_df$avg_log2FC.tf))
plot_data_df$avg_log2FC.tf.snRNA[plot_data_df$avg_log2FC.tf.snRNA == plot_data_df$TF] <- NA
plot_data_df$avg_log2FC.tf.snRNA <- as.numeric(plot_data_df$avg_log2FC.tf.snRNA)
summary(plot_data_df$diff)
summary(plot_data_df$x_plot)
## decide TFs to show
plot_text_right_df <- plot_data_df %>%
  filter(FDR < 0.05 & mean_score2 > 0 & !is.na(avg_log2FC.tf.snRNA) & avg_log2FC.tf.snRNA > 0.25 & !Is_FOSJUN & diff_2_vs_1 > 0)
tfs_right <- plot_text_right_df$TF
tfs_right <- c(tfs_right, "JUN", "TWIST1")
tfs_right

plot_text_left_df <- plot_data_df %>%
  filter(FDR < 0.05 & mean_score1 > 0 & !is.na(avg_log2FC.tf.snRNA) & avg_log2FC.tf.snRNA < -0.25 & diff_2_vs_1 < 0)
tfs_left <- plot_text_left_df$TF
# tfs_left <- tail(tfs_left, 25); tfs_left
## cap y axis
y_cap <- max(plot_data_df$Log10p_val_adj[!is.infinite(plot_data_df$Log10p_val_adj)])

# label TFS to show
plot_data_df <- plot_data_df %>%
  mutate(y_plot = ifelse(Log10p_val_adj >= y_cap, y_cap, Log10p_val_adj)) %>%
  mutate(TF_modified = gsub(x = TF, pattern = "\\(var.2\\)|\\(var.3\\)", replacement = "*")) %>%
  mutate(text_TF = ifelse(TF %in% c(tfs_right, tfs_left), TF_modified, NA))

# plot all markers--------------------------------------------------------------------
## plot
p <- ggplot()
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey70")
p <- p + geom_point_rast(data = subset(plot_data_df, y_plot < y_bottom), mapping = aes(x = x_plot, y = y_plot), alpha = 1, size = 1, color = "grey70", shape = 16)
p <- p + geom_point_rast(data = subset(plot_data_df, y_plot >= y_bottom & x_plot > 0), mapping = aes(x = x_plot, y = y_plot), alpha = 1, size = 1, color = color_right_deep, shape = 16)
p <- p + geom_point_rast(data = subset(plot_data_df, y_plot >= y_bottom & x_plot < 0), mapping = aes(x = x_plot, y = y_plot), alpha = 1, size = 1, color = color_left_deep, shape = 16)
# p <- p + scale_color_manual(values = c("FDR<0.05 (up)" = "red", "FDR<0.05 (down)" = "blue", "FDR<0.05" = "black", "FDR>=0.05" = "grey80"))
# p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_TF)),
#                          mapping = aes(x = x_plot, y = y_plot, label = text_TF), color = "black", force = 4, fontface = "bold", segment.alpha = 0.5)
p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_TF) & x_plot > 0),
                         mapping = aes(x = x_plot, y = y_plot, label = text_TF), max.overlaps = Inf,
                         color = "black", force = 4, alpha = 0.8, size = 5, segment.size = 0.2, segment.alpha = 0.5, xlim = c(0, 1.48))
p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_TF) & x_plot < 0), max.overlaps = Inf,
                         mapping = aes(x = x_plot, y = y_plot, label = text_TF),
                         color = "black", force = 5, alpha = 0.8, size = 5, segment.size = 0.2, segment.alpha = 0.5, xlim = c(-1.48, 0))
p <- p + ylim(0, 350)
p <- p + xlim(-1.5, 1.5)
p <- p + theme_classic()
p <- p + xlab("Motif score difference")
p <- p + ylab("-Log10(P-value-adjusted)")
p <- p + theme(axis.text = element_text(size = 14, color = "black"),
               axis.title = element_text(size = 15))
p
file2write <- paste0(dir_out, "volcano.", "png")
png(file2write, width = 800, height = 600, res = 150)
print(p)
dev.off()

file2write <- paste0(dir_out, "volcano.", "pdf")
pdf(file2write, width = 6, height = 5, useDingbats = F)
print(p)
dev.off()

