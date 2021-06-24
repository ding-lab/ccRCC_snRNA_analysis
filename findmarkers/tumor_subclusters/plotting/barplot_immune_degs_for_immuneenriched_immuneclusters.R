# Yige Wu @WashU Apr 2021

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
count_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/filter_degs/count_immune_degs_for_immuneenriched_tumorclusters/20210426.v1/Count.Immune_DEGs.tsv")

# prepare data ------------------------------------------------------------
plot_data_df <- count_df %>%
  filter(count_clusters >= 16/3)
plot_data_df$gene <- factor(x = plot_data_df$gene, levels = plot_data_df$gene)
# plot --------------------------------------------------------------------
p <- ggplot()
p <- p + geom_bar(data = plot_data_df, mapping = aes(x = gene, y = count_clusters, fill = mean_avg_logFC), stat =  "identity")
p <- p + theme_classic(base_size = 15)
p <- p + theme(axis.text.x = element_text(angle = 90, hjust=0.95,vjust=0.2))
p <- p + scale_fill_gradient(low = "yellow", high = "red", na.value = NA)
p <- p + theme(axis.title.x = element_blank())
file2write <- paste0(dir_out, "Count_clusters_with_DEGs.pdf")
pdf(file2write, width = 8, height = 3, useDingbats = F)
print(p)
dev.off()
