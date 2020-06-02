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
## input known CNV genes
knowncnvgenes_df <- readxl::read_xlsx(path = "./Resources/Known_Genetic_Alterations/Known_CNV.20200505.v1.xlsx", sheet = "Genes")

# make plot data  -----------------------------------------------------
## filter markers
plot_data_df <- markers_df %>%
  filter(p_val_adj < 0.05)
## add all clusters
plot_data_df <- merge(plot_data_df, 
                      tumorsubcluster_list_df %>%
                        filter(!is.na(Cluster_Manual)) %>%
                        select(Aliquot, Cluster_Manual),
                      by.x = c("id_aliquot", "cluster"),
                      by.y = c("Aliquot", "Cluster_Manual"), all = T)
plot_data_df <- unique(plot_data_df)
## add cluster name 
plot_data_df$id_aliquot_wu <- mapvalues(x = plot_data_df$id_aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
plot_data_df <- plot_data_df %>%
  mutate(deg_direction = ifelse(!is.na(avg_logFC),ifelse(avg_logFC > 0, "up", "down"), "NA")) %>%
  mutate(id_aliquot_cluster = paste0(id_aliquot_wu, "_C", (cluster + 1))) %>%
  mutate(name_tumorsubcluster = str_split_fixed(string = id_aliquot_cluster, pattern = "_", n = 2)[,2]) %>%
  mutate(x_plot = id_aliquot_cluster) %>%
  mutate(y_plot = avg_logFC)

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
## get top values and mark them in the data frame
plot_data_df$id_row <- seq.int(nrow(plot_data_df))
plot_text_df1 <- plot_data_df %>%
  group_by(id_aliquot_cluster) %>%
  top_n(1, avg_logFC)
plot_text_df2 <- plot_data_df %>%
  group_by(id_aliquot_cluster) %>%
  top_n(1, -avg_logFC)
plot_data_df$is_top <- (plot_data_df$id_row %in% c(as.vector(plot_text_df1$id_row), as.vector(plot_text_df2$id_row)))

# Create a jitter object for reproducible jitter:
jitter <- position_jitter(width = 0.2, height = 0.1)

# make plot ---------------------------------------------------------------
p <- ggplot()
p <- p + geom_point(data = plot_data_df, mapping = aes(x = x_plot, y = y_plot, size = (-log10(p_val_adj))^(1/3), fill = deg_direction), 
                    shape = 21, alpha = 0.3, position = jitter)
p <- p + facet_grid(.~id_aliquot_wu,
                    scales = "free", space = "free", shrink = T)
p <- p + scale_fill_manual(values = c("up" = "red", "down" = "blue", "NA" = "grey50"))
p <- p + geom_text_repel(data = plot_data_df, mapping = aes(x = x_plot, y = y_plot, label = ifelse(is_top, gene, NA)), angle = 90, size = 2.5)
p <- p + scale_x_discrete(breaks=plot_data_df$id_aliquot_cluster,
                          labels=plot_data_df$name_tumorsubcluster)
p <- p + ylim(c(-3,3))
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

