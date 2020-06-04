# Yige Wu @WashU May 2020

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

# input dependencies -------------------------------------------------------
## input DEGs
markers_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/findallmarkers_wilcox_tumorcells_by_manualsubcluster/20200427.v1/Tumormanualsubcluster.FindAllMarkers.Wilcox.Minpct0.1.Logfc0.25.tsv")
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200505.v1/meta_data.20200505.v1.tsv")
## input the complete list of subclusters
tumorsubcluster_list_df <- fread(data.table = F, input = "./Resources/snRNA_Processed_Data/Cell_Type_Assignment/Tumor_Subcluster/ccRCC_snRNA_Downstream_Processing - Individual.TumorCluster2Cell_Type.20200223.v1.tsv")
## input the bubble plot with the highlighted genes
cnv_plot <- readRDS(file = "./Resources/Analysis_Results/copy_number/summarize_cnv_fraction/bubble_cnv_fraction_per_tumormanualsubcluster/20200602.v2/bubbleplot_cnv_fraction.20200602.v2.RDS")
## input the cnv gene interactome
genes_interact_cnvgenes_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/filter_tumormanualsubcluster_deg_by_ppi/20200603.v1/tumormanualsubcluster.cnv_genes_interactome.20200603.v1.tsv")
## set p_val_adj_plot
p_val_adj_plot <- 0.05

# make plot data  -----------------------------------------------------
markers_df$id_aliquot_wu <- mapvalues(x = markers_df$id_aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
## filter markers
plot_data_df <- markers_df %>%
  filter(p_val_adj < p_val_adj_plot) %>%
  select(-id_aliquot_wu)
## add all clusters
plot_data_df <- merge(plot_data_df, 
                      tumorsubcluster_list_df %>%
                        filter(!is.na(Cluster_Manual)) %>%
                        select(Aliquot, Cluster_Manual),
                      by.x = c("id_aliquot", "cluster"),
                      by.y = c("Aliquot", "Cluster_Manual"), all = T)
plot_data_df <- unique(plot_data_df)
plot_data_df$id_aliquot_wu <- mapvalues(x = plot_data_df$id_aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
## add cluster name 
plot_data_df <- plot_data_df %>%
  mutate(deg_direction = ifelse(!is.na(avg_logFC), 
                                ifelse(avg_logFC > 0, 
                                       ifelse(p_val_adj < p_val_adj_plot, paste0("logFC>0 & p_val_adj<", p_val_adj_plot), "logFC>0 & p_val<0.05"), 
                                       ifelse(p_val_adj < p_val_adj_plot, paste0("logFC<0 & p_val_adj<", p_val_adj_plot), "logFC<0 & p_val < 0.05")), 
                                "NA")) %>%
  mutate(id_aliquot_cluster = paste0(id_aliquot_wu, "_C", (cluster + 1))) %>%
  mutate(name_tumorsubcluster = str_split_fixed(string = id_aliquot_cluster, pattern = "_", n = 2)[,2]) %>%
  mutate(x_plot = id_aliquot_cluster) %>%
  mutate(y_plot = avg_logFC)
plot_data_df$id_row <- seq.int(nrow(plot_data_df))

# make plot dependencies --------------------------------------------------
## count the number of subclusters per aliquot
count_subclusters_df <- plot_data_df %>%
  select(id_aliquot_wu, id_aliquot_cluster) %>%
  unique() %>%
  select(id_aliquot_wu) %>%
  table() %>%
  as.data.frame() %>%
  rename(id_aliquot_wu = '.') %>%
  rename(count_subclusters = Freq) %>%
  arrange(count_subclusters)
levels_id_aliquot_wu <- count_subclusters_df$id_aliquot_wu
## order aliquot id by transformating into factor
plot_data_df$id_aliquot_wu <- factor(plot_data_df$id_aliquot_wu, levels = levels_id_aliquot_wu)
## show the genes highlighted in the cnv plot
cnv_plot_data_df <- cnv_plot$data
cnv_hl_data_df1 <- cnv_plot_data_df %>%
  filter(!is.na(x2)) %>%
  select(aliquot.wu, gene_symbol) %>%
  unique() %>%
  mutate(is_highlight = T)
plot_data_df <- merge(plot_data_df, cnv_hl_data_df1, by.x = c("id_aliquot_wu", "gene"), by.y = c("aliquot.wu", "gene_symbol"), all.x = T)
genes_interact_cnvgenes_top_df <- genes_interact_cnvgenes_df %>%
  group_by(name_cluster, cna_gene_symbol) %>%
  top_n(3, avg_logFC)
cnv_hl_data_df2 <- genes_interact_cnvgenes_top_df %>%
  as.data.frame() %>%
  rename(id_aliquot_cluster = name_cluster) %>%
  rename(aliquot.wu = id_aliquot_wu) %>%
  rename(gene_symbol = gene) %>%
  select(id_aliquot_cluster, gene_symbol) %>%
  unique() %>%
  mutate(is_highlight_cnvgeneinteractome = T)
plot_data_df <- merge(plot_data_df, cnv_hl_data_df2, by.x = c("id_aliquot_cluster", "gene"), by.y = c("id_aliquot_cluster", "gene_symbol"), all.x = T)
plot_data_df$is_highlight[is.na(plot_data_df$is_highlight)] <- F
plot_data_df$is_highlight_cnvgeneinteractome[is.na(plot_data_df$is_highlight_cnvgeneinteractome)] <- F
# Create a jitter object for reproducible jitter:
jitter <- position_jitter(width = 0.2, height = 0.1)
## make colors for different levels of DEGs
colors_deg_direction <- RColorBrewer::brewer.pal(n = 5, name = "PuOr")
names(colors_deg_direction) <- c(paste0("logFC>0 & p_val_adj<", p_val_adj_plot), "logFC>0 & p_val<0.05",
                                 "NA",
                                 "logFC<0 & p_val < 0.05", paste0("logFC<0 & p_val_adj<", p_val_adj_plot))
# make plot ---------------------------------------------------------------
p <- ggplot()
# p <- p + geom_point(data = plot_data_df, mapping = aes(x = x_plot, y = y_plot, size = (-log10(p_val))^(1/3), fill = deg_direction), 
p <- p + geom_point(data = plot_data_df, mapping = aes(x = x_plot, y = y_plot, size = (-log10(p_val_adj))^(1/3), fill = deg_direction),
                    shape = 21, alpha = 0.3, position = jitter)
p <- p + facet_grid(.~id_aliquot_wu,
                    scales = "free", space = "free", shrink = T)
p <- p + scale_fill_manual(values = colors_deg_direction)
# p <- p + geom_text_repel(data = plot_data_df, mapping = aes(x = x_plot, y = y_plot, label = ifelse(is_highlight, gene, NA)), angle = 90, size = 3, color = "black", ylim = c(1, NA))
# p <- p + geom_text_repel(data = subset(plot_data_df, avg_logFC > 0), mapping = aes(x = x_plot, y = y_plot, label = ifelse(is_highlight | is_highlight_cnvgeneinteractome, gene, NA)), angle = 90, size = 3, color = "black", ylim = c(1, NA))
p <- p + scale_x_discrete(breaks=plot_data_df$id_aliquot_cluster,
                          labels=plot_data_df$name_tumorsubcluster)
p <- p + ylim(c(-2.25,2.25))
p <- p + theme_bw()
p <- p + theme(panel.spacing = unit(0, "lines"))
p <- p + theme(strip.text.x = element_text(angle = 90))
p
# save output -------------------------------------------------------------
## save plot as an object
file2write <- paste0(dir_out, "bubbleplot_deg_by_tumorsubcluster.", run_id, ".RDS")
saveRDS(object = p, file = file2write, compress = T)
## save plot
png(file = paste0(dir_out, "bubbleplot_deg_by_tumorsubcluster.", run_id, ".png"), 
    width = 2500, height = 1000, res = 150)
print(p)
dev.off()

