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
dam_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Enriched_Motifs/Score_difference.Tumor_Normal_comparison.20210509.tsv")

# set plotting parameters -------------------------------------------------
## set y bottom threshold
y_bottom <- -log10(0.05)
## set x limits to distinguish colors
x_pos <- 0.5
x_neg <- -0.5
## colors
color_purple <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[4]
color_red <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[1]
color_blue <- RColorBrewer::brewer.pal(n = 4, name = "Set1")[2]

# make data for plotting --------------------------------------------------
dam_sig_df <- dam_df %>%
  mutate(log10_pvalue = -log10(pvalue)) %>%
  mutate(log10_pvalue_capped = ifelse(is.infinite(log10_pvalue), 150, log10_pvalue)) %>%
  filter(FDR < 0.05)

plot_data_df <- dam_sig_df %>%
  group_by(TF_Name) %>%
  summarise(Num_sig_up = length(which(diff > 0)), Num_sig_down = length(which(diff < 0)),
            avg_log10_pvalue = mean(log10_pvalue_capped), avg_diff = mean(diff)) %>%
  mutate(foldchange_type = ifelse(Num_sig_down == 0, "consistently higher in ccRCC",
                                  ifelse(Num_sig_up == 0, "consistently lower in ccRCC", "mixed fold change directions"))) %>%
  # mutate(number_foldchange = ifelse(avg_diff > 0, Num_sig_up, Num_sig_down)) %>%
  mutate(size_plot = abs(Num_sig_up - Num_sig_down)/24) %>%
  mutate(x_plot = ifelse(avg_diff < -2, -2, ifelse(avg_diff > 2, 2,  avg_diff))) %>%
  mutate(y_plot = avg_log10_pvalue) %>%
  arrange(desc(foldchange_type))
## decide TFs to show
tfnames_show <- plot_data_df$TF_Name[plot_data_df$size_plot == 1]
## label TFS to show
plot_data_df <- plot_data_df %>%
  # mutate(TF_modified = gsub(x = TF_Name, pattern = "\\(var.2\\)", replacement = "*")) %>%
  mutate(text_tf = ifelse(TF_Name %in% tfnames_show, TF_Name, NA))

# plot all markers--------------------------------------------------------------------
## plot
p <- ggplot()
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey70")
p <- p + geom_point(data = subset(plot_data_df, foldchange_type == "mixed fold change directions"), mapping = aes(x = x_plot, y = y_plot, size = size_plot, color = foldchange_type), alpha = 0.5)
p <- p + geom_point(data = subset(plot_data_df, foldchange_type != "mixed fold change directions"), mapping = aes(x = x_plot, y = y_plot, size = size_plot, color = foldchange_type), alpha = 0.9)
# p <- p + geom_point(data = subset(plot_data_df, foldchange_type == "Mixed fold change directions"), mapping = aes(x = x_plot, y = y_plot, size = size_plot), alpha = 0.5, color = color_purple)
# p <- p + geom_point(data = subset(plot_data_df, foldchange_type == "Consistently higher in ccRCC"), mapping = aes(x = x_plot, y = y_plot, size = size_plot), alpha = 0.8, color = color_red)
# p <- p + geom_point(data = subset(plot_data_df, foldchange_type == "Consistently lower in ccRCC"), mapping = aes(x = x_plot, y = y_plot, size = size_plot), alpha = 0.8, color = color_blue)
p <- p + scale_color_manual(values = c("consistently higher in ccRCC" = color_red, 
                                       "consistently lower in ccRCC" = color_blue, 
                                       "mixed fold change directions" = color_purple))
# p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_TF)),
#                          mapping = aes(x = x_plot, y = y_plot, label = text_TF), color = "black", force = 4, fontface = "bold", segment.alpha = 0.5)
p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_tf) & x_plot > 0),
                         mapping = aes(x = x_plot, y = y_plot, label = text_tf), 
                         color = "black", alpha = 1, size = 3, fontface = "bold",
                         segment.size = 0.4, segment.alpha = 1, min.segment.length = 0,
                         xlim = c(0, NA))
p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_tf) & x_plot < 0),
                         mapping = aes(x = x_plot, y = y_plot, label = text_tf), 
                         color = "black", alpha = 1, size = 3, fontface = "bold",
                         segment.size = 0.4, segment.alpha = 1, min.segment.length = 0,
                         xlim = c(NA, 0))
p <- p + scale_size_area(max_size = 4)
p <- p + theme_classic()
p <- p + xlab("Motif score difference (ccRCC cells - PT cells)")
p <- p + ylab("-Log10FDR")
p <- p + guides(color = guide_legend(title = "Motif type", title.position = "top", nrow = 3, override.aes = aes(size = 3)),
                size = guide_legend(title = "Motif score difference\nconsistency index", title.position = "top", nrow = 2))
# p <- p + labs(color = "|(No. tumors with higher motif scores) - (No. tumors with lower motif scores)|/(No. all tumors)")
p <- p + theme(axis.text = element_text(size = 14),
               axis.title = element_text(size = 15),
               legend.position = "bottom", legend.box = "horizontal")
file2write <- paste0(dir_out, "volcano.", "png")
png(file2write, width = 800, height = 800, res = 150)
print(p)
dev.off()
# file2write <- paste0(dir_out, "volcano.", "pdf")
# pdf(file2write, width = 6, height = 5, useDingbats = F)
# print(p)
# dev.off()