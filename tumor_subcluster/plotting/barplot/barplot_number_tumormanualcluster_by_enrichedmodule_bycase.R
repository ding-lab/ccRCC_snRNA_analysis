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
## input meta data
id_metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210322.v1/meta_data.20210322.v1.tsv")
## input by cluster enrichment assignment
enrich_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/assign_tumorcluster_by_msigdb_geneset_scores/20210421.v3/MsigDB_Hallmark.Top15GeneSets.4Module.Enrichment.tsv")

# make plot data ----------------------------------------------------------
plotdata_df <- enrich_df %>%
  mutate(sample = str_split_fixed(string = cluster_name, pattern = "_", n = 2)[,1]) %>%
  mutate(case = paste0(str_split_fixed(string = sample, pattern = "\\.", n = 3)[,1], "-", str_split_fixed(string = sample, pattern = "\\.", n = 3)[,2])) %>%
  mutate(cluster_enrich_type = ifelse(EMT, "EMT",
                                      ifelse(Cell_cycle, "Cell_cycle",
                                      ifelse(Immune, "Immune",
                                             ifelse(mTOR, "mTOR", "Other")))))
x_anno_df <- plotdata_df %>%
  select(sample) %>%
  table() %>%
  as.data.frame() %>%
  rename(sample = '.')
plotdata_df$x_plot_order <- mapvalues(x = plotdata_df$sample, from = x_anno_df$sample, to = as.vector(x_anno_df$Freq))
plotdata_df$x_plot_order <- as.numeric(plotdata_df$x_plot_order)
plotdata_df <- plotdata_df %>%
  arrange(desc(x_plot_order), cluster_enrich_type)
plotdata_df$x_plot <- factor(x = plotdata_df$sample, levels = unique(plotdata_df$sample))
## make colors
colors_enrich_type <- RColorBrewer::brewer.pal(n = 8, name = "Set2")[c(1,6,3,4,8)]
names(colors_enrich_type) <- c("EMT", "Cell_cycle", "Immune", "mTOR", "Other")
# plot --------------------------------------------------------------------
p <- ggplot()
p <- p + geom_bar(data = plotdata_df, mapping = aes(x = x_plot, y = 1, fill = cluster_enrich_type), stat = "identity")
p <- p + scale_fill_manual(values = colors_enrich_type)
p <- p + theme_classic(base_size = 15)
p <- p + theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.2))
p <- p + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank())
p <- p + scale_y_continuous(breaks = seq(0, 10, 2))
p <- p + ylab("Number of tumor subclusters / sample")
file2write <- paste0(dir_out, "barplot.pdf")
pdf(file2write, width = 8, height = 5, useDingbats = F)
print(p)
dev.off()



