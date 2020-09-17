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
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/findallmarker_wilcox_tumor_vs_pt_on_katmai/20200903.v1/findallmarkers_wilcox_tumorcells_vs_pt.20200903.v1.tsv")
## input the peaks to TF motifs
# peaks2motif_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/annotate_peaks/annotate_tumor_vs_normal_da_peaks_to_celltypespecific_tfmotifs/20200916.v1/DA_peaks.Filter_for_Tumor_vs_Normal.MotifAnnotated.20200916.v1.tsv")
deg2tfgene_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/annotate_deg/annotate_tumor_vs_normal_degs_to_tf_by_snatac/20200916.v1/DEGs_with_DA_peaks.20200916.v1.tsv")

# for (motifname_tmp in c("HIF1A")) {
for (motifname_tmp in unique(peaks2motif_df$motif.name)) {
  ## get the cell type and rank for the motif
  celltype_tmp <- peaks2motif_df$Cell_type.filename[peaks2motif_df$motif.name == motifname_tmp] %>% unique()
  rank_motif_tmp <- peaks2motif_df$rank_motif_by_avglogFC[peaks2motif_df$motif.name == motifname_tmp] %>% unique()
  
  ## set output directory
  dir_out_tmp <- paste0(dir_out, celltype_tmp, "/")
  dir.create(dir_out_tmp)
  
  # set genes to filter to show -----------------------------------------------------
  genes_filter <- unique(peaks2motif_df$SYMBOL[peaks2motif_df$motif.name == motifname_tmp])
  
  # make data for plotting --------------------------------------------------
  plot_data_df <- deg_df %>%
    dplyr::rename(genesymbol = row_name) %>%
    mutate(Log10p_val_adj = -log10(x = p_val_adj)) %>%
    mutate(x_plot = avg_logFC) %>%
    mutate(text_gene = ifelse(genesymbol %in% genes_filter, genesymbol, NA))
  ## cap y axis
  y_cap <- max(plot_data_df$Log10p_val_adj[!is.infinite(plot_data_df$Log10p_val_adj)])
  plot_data_df <- plot_data_df %>%
    mutate(y_plot = ifelse(Log10p_val_adj >= y_cap, y_cap, Log10p_val_adj))
  ## set y top shown
  y_top_area <- y_cap*1.1
  ## set y bottom threshold
  y_bottom <- -log10(0.05)
  ## set x limits to distinguish colors
  x_pos <- log(2)
  x_neg <- -log(2)
  
  # plot all markers--------------------------------------------------------------------
  ## plot
  p <- ggplot()
  # p <- p + geom_hline(yintercept = y_bottom, linetype = 2)
  p <- p + geom_point(data = subset(plot_data_df, x_plot >= x_pos), mapping = aes(x = x_plot, y = y_plot), alpha = 0.8, size = 0.5, color = "red")
  p <- p + geom_point(data = subset(plot_data_df, x_plot <= x_neg), mapping = aes(x = x_plot, y = y_plot), alpha = 0.8, size = 0.5, color = "blue")
  # p <- p + scale_color_manual(values = c("FDR<0.05 (up)" = "red", "FDR<0.05 (down)" = "blue", "FDR<0.05" = "black", "FDR>=0.05" = "grey80"))
  p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_gene)),
                           mapping = aes(x = x_plot, y = y_plot, label = text_gene), color = "black", force = 4, fontface = "bold", segment.alpha = 0.5)
  # p <- p + geom_text_repel(data = subset(plot_data_df, !is.na(text_gene) & avg_logFC <= 0),
  #                          mapping = aes(x = avg_logFC, y = y_capped, label = text_gene), color = "black", force = 2)
  p <- p + theme_bw()
  p <- p + ggtitle(subtitle = paste0("And differentially expressed between tumor vs normal pt cells (snRNA)"), 
                   label = paste0("Genes with ", motifname_tmp, " motif in the promoter region differentially\naccessible between tumor vs normal pt cells (snATAC)"))
  p <- p + xlim(c(-3, 3))
  p <- p + ylim(c(0, y_top_area))
  p <- p + xlab("ln(Fold-Change) (tumor cells vs normal proximal tubule cells)")
  p <- p + ylab("-Log10(P-value-adjusted)")
  p
  file2write <- paste0(dir_out_tmp, celltype_tmp, "_motif_rank", rank_motif_tmp, ".",  motifname_tmp, ".png")
  png(file2write, width = 800, height = 800, res = 150)
  print(p)
  dev.off()
}
