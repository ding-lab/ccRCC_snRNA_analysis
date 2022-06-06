# Yige Wu @WashU May 2022

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "ggplot2"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    print(paste0(pkg_name_tmp, "is being installed!"))
    BiocManager::install(pkgs = pkg_name_tmp, update = F)
    install.packages(pkg_name_tmp, dependencies = T)
  }
  print(paste0(pkg_name_tmp, " is installed!"))
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
exp_df <- fread(input = "./Resources/Analysis_Results/bulk/expression/rna/knockout_cell_lines/normalization/get_edgeR_TMM_normalized_counts_all_cell_lines/20220603.v1/CPM.TMM_normalized.All_Cell_Lines.20220603.v1.tsv", data.table = F)

# specify parameters ------------------------------------------------------
colnames_value <- colnames(exp_df)[grepl(pattern = "sample", x = colnames(exp_df))]
exp_df <- as.data.table(exp_df)
# genes_plot <- c("KLF9", "CP")

genes_plot <- c("UCHL1")
# genes_plot <- c("VEGFA")

# format expression data --------------------------------------------------
plot_data_long_df <- exp_df %>%
  filter(external_gene_name %in% genes_plot) %>%
  melt.data.table(measure.vars = colnames_value) %>%
  mutate(sample = gsub(x = variable, pattern = "sample\\.", replacement = "")) %>%
  mutate(sample_text = sample)

p <- ggplot()
p <- p + geom_col(data = plot_data_long_df, mapping = aes(x = sample_text, y = value), position=position_dodge())
p <- p + theme_classic()
p <- p + ylab(label = "CPM")
p <- p + ggtitle(label = paste0(paste0(genes_plot, collapse = " & "), " expression")) # , subtitle = "by RNA-seq"
# p <- p + theme(strip.background = element_rect(fill = NA),
#                panel.spacing = unit(0, "lines"))
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p <- p + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank())
p
file2write <- paste0(dir_out, paste0(genes_plot, collapse = "_"), ".bulkRNA.CPM.", "pdf")
pdf(file2write, width = 2, height = 2.5, useDingbats = F)
print(p)
dev.off()
file2write <- paste0(dir_out, paste0(genes_plot, collapse = "_"), ".bulkRNA.CPM.", "png")
png(file2write, width = 300, height = 450, res = 150)
print(p)
dev.off()

## when genes_plot == "CP
exp_shRNA <- sum(plot_data_long_df$value[plot_data_long_df$sample %in% c("caki_1_cp_c2_e1", "caki_1_cp_c1_e1")])
exp_NT <- sum(plot_data_long_df$value[plot_data_long_df$sample %in% c("caki_1_control_e1", "dr_caki_1_rna")])
exp_shRNA/exp_NT

