# Yige Wu @WashU Jun 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "data.table",
  "stringr",
  "ggpubr",
  "ggplot2"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    install.packages(pkg_name_tmp, dependencies = T)
  }
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg_name_tmp)
  }
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the tumor clusters
clusters_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/count/count_cellnumber_per_manual_cluster_rm_doublet/20220307.v1/CellNumberPerTumorManualCluster.20220307.v1.tsv")

# plot by sample --------------------------------------------------------------------
plot_data_df <- clusters_df %>%
  mutate(case = str_split_fixed(string = easy_id, pattern = "\\-T", n = 2)[,1]) %>%
  filter(Freq >= 50) %>%
  filter(case != "C3L-00359") %>%
  group_by(easy_id) %>%
  summarise(Freq = n())
p <- ggplot(data = plot_data_df, mapping = aes(x = Freq))
p <- p + geom_histogram(color="black", fill="white", binwidth = 1)
p <- p + scale_x_continuous(breaks = 1:10)
p <- p + xlab("Number of Tumor Subclusters")
p <- p + ylab("Number of Samples")
p <- p + theme_bw()
p <- p + theme(axis.text.y = element_text(size = 15, color = "black"),
               axis.text.x = element_text(vjust = 0.1, size = 15, color = "black"))
p
pdf2write <- paste0(dir_out, "NumberOfManualTumorSubclusters.persample.", ".pdf")
pdf(file = pdf2write, width = 4, height = 3, useDingbats = F)
print(p)
dev.off()


# plot by sample --------------------------------------------------------------------
plot_data_df <- clusters_df %>%
  mutate(case = str_split_fixed(string = easy_id, pattern = "\\-T", n = 2)[,1]) %>%
  filter(Freq >= 50) %>%
  filter(case != "C3L-00359") %>%
  group_by(case) %>%
  summarise(Freq = n())
p <- ggplot(data = plot_data_df, mapping = aes(x = Freq))
p <- p + geom_histogram(color="black", fill="white", binwidth = 1)
p <- p + scale_x_continuous(breaks = 1:10)
p <- p + xlab("Number of Tumor Subclusters")
p <- p + ylab("Number of Patients")
p <- p + theme_bw()
p <- p + theme(axis.text.y = element_text(size = 15, color = "black"),
               axis.text.x = element_text(vjust = 0.1, size = 15, color = "black"))
p
pdf2write <- paste0(dir_out, "NumberOfManualTumorSubclusters.perpatient.", ".pdf")
pdf(file = pdf2write, width = 8, height = 5, useDingbats = F)
print(p)
dev.off()