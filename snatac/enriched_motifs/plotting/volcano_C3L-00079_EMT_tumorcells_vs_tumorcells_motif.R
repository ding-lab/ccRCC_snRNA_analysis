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
## input degs
da_df <- fread(data.table = F, input = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/Enriched_Motifs/TransitionalPopulation/C3L_00079_Transitional_TumorCells_Score_difference.Wilcoxon.two_sided.20201113.tsv")
## 
da_bycelltype_df <- fread(data.table = F, input = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/Enriched_Motifs/TransitionalPopulation/C3L_00079_T1_Score_difference.20201111.tsv")

# filter the emt genes ----------------------------------------------------
da_bycelltype_filtered_df <- da_bycelltype_df %>%
  filter(p_val_adj < 0.05) %>%
  arrange(desc(diff))
da_bycelltype_filtered_df <- da_bycelltype_filtered_df[!duplicated(da_bycelltype_filtered_df$TF),]
tfs_emt_tumorcells <- da_bycelltype_filtered_df$TF[da_bycelltype_filtered_df$cell_t2 == "Transitional cells"]
tfs_tumorcells <- da_bycelltype_filtered_df$TF[da_bycelltype_filtered_df$cell_t2 == "Tumor cells"]

# set plotting parameters -------------------------------------------------
## set y bottom threshold
y_bottom <- -log10(0.05)
## set x limits to distinguish colors
x_pos <- 0.5
x_neg <- -0.5
## colors
color_right_deep <- RColorBrewer::brewer.pal(n = 12, name = "Paired")[12]
color_right_pale <- RColorBrewer::brewer.pal(n = 7, name = "Set2")[7]
color_left_deep <- RColorBrewer::brewer.pal(n = 6, name = "Dark2")[4]
color_left_pale <- RColorBrewer::brewer.pal(n = 6, name = "Set2")[4]

# make data for plotting --------------------------------------------------
plot_data_df <- da_df %>%
  mutate(Is_FOSJUN = grepl(pattern = "FOS|JUN", x = TF)) %>%
  mutate(Log10p_val_adj = -log10(x = p_val_adj)) %>%
  mutate(x_plot = diff) %>%
  arrange(desc(diff))
## cap y axis
y_cap <- max(plot_data_df$Log10p_val_adj[!is.infinite(plot_data_df$Log10p_val_adj)])
## decide TFs to show
plot_data_filtered_df <- plot_data_df %>%
  filter(TF %in% c(tfs_emt_tumorcells, tfs_tumorcells))
tfs_right <- plot_data_filtered_df$TF[!(plot_data_filtered_df$Is_FOSJUN) & plot_data_filtered_df$diff > 0 & plot_data_filtered_df$mean_score1 > 0 & plot_data_filtered_df$TF %in% tfs_emt_tumorcells] %>%
  head(19)
tfs_right <- c(tfs_right, "FOS::JUN")
tfs_left <- plot_data_filtered_df$TF[!(plot_data_filtered_df$Is_FOSJUN) & plot_data_filtered_df$diff < 0 & plot_data_filtered_df$mean_score2 > 0 & plot_data_filtered_df$TF %in% tfs_tumorcells] %>%
  tail(20)
tfs_left
## label TFS to show
plot_data_df <- plot_data_df %>%
  mutate(y_plot = ifelse(Log10p_val_adj >= y_cap, y_cap, Log10p_val_adj)) %>%
  mutate(TF_modified = gsub(x = TF, pattern = "\\(var.2\\)", replacement = "*")) %>%
  mutate(text_TF = ifelse(TF %in% c(tfs_right, tfs_left), TF_modified, NA))

# plot all markers--------------------------------------------------------------------
## plot
p <- ggplot()
# p <- p + geom_vline(xintercept = x_pos, linetype = 2)
# p <- p + geom_vline(xintercept = x_neg, linetype = 2)
p <- p + geom_vline(xintercept = 0, linetype = 2)
p <- p + geom_point(data = subset(plot_data_df, y_plot < y_bottom), mapping = aes(x = x_plot, y = y_plot), alpha = 0.5, size = 0.5, color = "grey70")
p <- p + geom_point(data = subset(plot_data_df, y_plot >= y_bottom & x_plot > 0), mapping = aes(x = x_plot, y = y_plot), alpha = 0.5, size = 0.5, color = color_right_deep)
p <- p + geom_point(data = subset(plot_data_df, y_plot >= y_bottom & x_plot < 0), mapping = aes(x = x_plot, y = y_plot), alpha = 0.5, size = 0.5, color = color_left_deep)
# p <- p + scale_color_manual(values = c("FDR<0.05 (up)" = "red", "FDR<0.05 (down)" = "blue", "FDR<0.05" = "black", "FDR>=0.05" = "grey80"))
# p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_TF)),
#                          mapping = aes(x = x_plot, y = y_plot, label = text_TF), color = "black", force = 4, fontface = "bold", segment.alpha = 0.5)
p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_TF) & x_plot > 0),
                         mapping = aes(x = x_plot, y = y_plot, label = text_TF), color = "black", force = 4, fontface = "italic", segment.alpha = 1, size = 6, segment.size = 0.2)
p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_TF) & x_plot < 0),
                         mapping = aes(x = x_plot, y = y_plot, label = text_TF), color = "black", force = 3, fontface = "italic", segment.alpha = 1, segment.size = 0.2)
p <- p + xlim(-3, NA)
p <- p + theme_bw()
# p <- p + ggtitle(label = paste0("C3L-00079 ", "VIM-high transitional cells vs tumor cells"))
p <- p + xlab("Bias-corrected accessibility deviations")
p <- p + ylab("-Log10(P-value-adjusted)")
p <- p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
               panel.background = element_blank())
# p <- p + theme(axis.line = element_line(colour = "black"))
p
file2write <- paste0(dir_out, "volcano.", "png")
png(file2write, width = 800, height = 600, res = 150)
print(p)
dev.off()

file2write <- paste0(dir_out, "volcano.", "pdf")
pdf(file2write, width = 6, height = 5, useDingbats = F)
print(p)
dev.off()