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

plot_data_df <- dam_df %>%
  mutate(log10_pvalue = -log10(pvalue)) %>%
  mutate(log10_pvalue_capped = ifelse(is.infinite(log10_pvalue), 150, log10_pvalue)) %>%
  group_by(TF_Name) %>%
  summarise(Num_sig_up = length(which(diff > 0 & FDR < 0.05)), Num_sig_down = length(which(diff < 0 & FDR < 0.05)),
            Num_up = length(which(diff > 0)), Num_down = length(which(diff < 0)),
            avg_log10_pvalue = mean(log10_pvalue_capped), avg_diff = mean(diff)) %>%
  mutate(foldchange_type = ifelse(Num_down == 0 & Num_sig_up >= 12, "consistently higher in ccRCC",
                                  ifelse(Num_up == 0 & Num_sig_down >= 12, "consistently lower in ccRCC", "insignificant"))) %>%
  mutate(x_plot = ifelse(avg_diff < -2, -2, ifelse(avg_diff > 2, 2,  avg_diff))) %>%
  mutate(y_plot = avg_log10_pvalue) %>%
  arrange(desc(foldchange_type))
## decide TFs to show
tfnames_show <- plot_data_df$TF_Name[plot_data_df$Num_sig_up == 24 | plot_data_df$Num_sig_down == 24]
## label TFS to show
plot_data_df <- plot_data_df %>%
  # mutate(TF_modified = gsub(x = TF_Name, pattern = "\\(var.2\\)", replacement = "*")) %>%
  mutate(text_tf = ifelse(TF_Name %in% tfnames_show, TF_Name, NA))

# plot all markers--------------------------------------------------------------------
p <- ggplot()
p <- p + geom_vline(xintercept = 0, linetype = 2, color = "grey70")
p <- p + geom_point(data = subset(plot_data_df, foldchange_type == "insignificant"), 
                    mapping = aes(x = x_plot, y = y_plot, color = foldchange_type), alpha = 0.5, shape = 16)
p <- p + geom_point(data = subset(plot_data_df, foldchange_type != "insignificant"), 
                    mapping = aes(x = x_plot, y = y_plot, color = foldchange_type), alpha = 0.9, shape = 16, size = 3)
p <- p + scale_color_manual(values = c("consistently higher in ccRCC" = color_red, 
                                       "consistently lower in ccRCC" = color_blue, 
                                       "insignificant" = "grey50"))
p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_tf) & x_plot > 0),
                         mapping = aes(x = x_plot, y = y_plot, label = text_tf), 
                         color = "black", alpha = 1, size = 5, #fontface = "bold",
                         segment.size = 0.2, segment.alpha = 0.5, min.segment.length = 0,
                         xlim = c(0, NA))
p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_tf) & x_plot < 0),
                         mapping = aes(x = x_plot, y = y_plot, label = text_tf), 
                         color = "black", alpha = 1, size = 5, #fontface = "bold",
                         segment.size = 0.2, segment.alpha = 0.5, min.segment.length = 0,
                         xlim = c(NA, 0), ylim = c(75, NA))
p <- p + scale_size_area(max_size = 4)
p <- p + theme_classic()
p <- p + xlab("Motif score difference for\nccRCC cells vs. PT cells")
p <- p + ylab("-Log10FDR")
p <- p + theme(axis.text = element_text(size = 14, color = "black"),
               axis.title = element_text(size = 14),
               legend.position = "none")

file2write <- paste0(dir_out, "volcano.nolegend.", "pdf")
pdf(file2write, width = 5.5, height = 5, useDingbats = F)
print(p)
dev.off()
