# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
# remotes::install_github("coolbutuseless/ggpattern")
library(ggpattern)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input by cluster enrichment assignment
enrich_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/assign_tumorcluster_by_msigdb_geneset_scores/20210805.v1/MsigDB_Hallmark.Top15GeneSets.4Module.Enrichment.tsv")

# make plot data ----------------------------------------------------------
plotdata_df <- enrich_df %>%
  mutate(sample = str_split_fixed(string = cluster_name, pattern = "_", n = 2)[,1]) %>%
  mutate(sample = gsub(x = sample, pattern = "\\.", replacement = "-")) %>%
  mutate(case = paste0(str_split_fixed(string = sample, pattern = "\\.", n = 3)[,1], "-", str_split_fixed(string = sample, pattern = "\\.", n = 3)[,2])) %>%
  mutate(cluster_enrich_type = ifelse(EMT, "EMT",
                                      ifelse(Cell_cycle, "Cell_cycle",
                                             ifelse(Immune, "Immune",
                                                    ifelse(mTOR, "mTOR", "Other"))))) %>%
  mutate(cluster_enrich_type2 = ifelse(cluster_enrich_type != "Cell_cycle" & Cell_cycle, "Cell_cycle",
                                       ifelse(cluster_enrich_type != "Immune" & Immune, "Immune",
                                              ifelse(cluster_enrich_type != "mTOR" & mTOR, "mTOR", "NA")))) %>%
  mutate(cluster_enrich_type3 = ifelse(cluster_enrich_type != "Immune" & cluster_enrich_type2 != "Immune" & Immune, "Immune",
                                       ifelse(cluster_enrich_type != "mTOR" & cluster_enrich_type2 != "mTOR" & mTOR, "mTOR", "NA")))
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
p <- ggplot(data = plotdata_df, mapping = aes(x_plot, 1))
p <- p + geom_col_pattern(
  aes(fill = cluster_enrich_type, pattern_fill = cluster_enrich_type2, pattern_density = cluster_enrich_type2, pattern_size = 0),
  colour          = 'black',
  pattern         = 'stripe'
)
p <- p + scale_fill_manual(values = colors_enrich_type)
p <- p + scale_pattern_fill_manual(values = colors_enrich_type)
p <- p + scale_pattern_density_manual(values = c("NA" = 0, "Cell_cycle"=0.5, "Immune"=0.5, "mTOR" = 0.5))
p <- p + theme_classic(base_size = 17)
p <- p + theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.2, size = 17), axis.line.x = element_blank(), axis.ticks.x = element_blank())
p <- p + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.line.y = element_blank())
p <- p + theme(axis.text.y = element_text(size = 17), axis.ticks.y = element_blank())
p <- p + scale_y_continuous(breaks = seq(0, 10, 2))
p <- p + ylab("Number of tumor subclusters/sample")
file2write <- paste0(dir_out, "barplot.pdf")
pdf(file2write, width = 10, height = 6, useDingbats = F)
print(p)
dev.off()


p <- ggplot(data = plotdata_df, mapping = aes(x_plot, 1))
p <- p + geom_col_pattern(
  aes(fill = cluster_enrich_type, pattern_fill = cluster_enrich_type2, pattern_density = cluster_enrich_type2),
  colour          = 'black',
  pattern         = 'stripe'
)
p <- p + scale_fill_manual(values = colors_enrich_type)
p <- p + scale_pattern_fill_manual(values = colors_enrich_type)
p <- p + scale_pattern_density_manual(values = c("NA" = 0, "Cell_cycle"=0.5, "Immune"=0.5, "mTOR" = 0.5))
p <- p + theme_classic(base_size = 15)
p <- p + theme(axis.text.x = element_text(angle = 90,hjust=0.95,vjust=0.2))
p <- p + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank())
p <- p + scale_y_continuous(breaks = seq(0, 10, 2))
p <- p + ylab("Number of tumor subclusters\nper sample")
file2write <- paste0(dir_out, "barplot.png")
png(file2write, width = 1100, height = 600, res = 150)
print(p)
dev.off()

