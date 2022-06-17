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
exp_df <- fread(input = "./Resources/Analysis_Results/bulk/expression/rna/knockout_cell_lines/divide_knockout_cellline_TPM_TMMnormalized_byscrambled/20220516.v1/Knockout_celllines.Kallisto.gene_level.CPM.TMMnormalized.DividedByControl.20220516.v1.tsv", data.table = F)

# specify parameters ------------------------------------------------------
colnames_value <- colnames(exp_df)[grepl(pattern = "sample", x = colnames(exp_df))]
exp_df <- as.data.table(exp_df)
# genes_plot <- c("KLF9", "CP")
genes_plot <- c("MXI1", "CP")
genes_plot <- c("MXI1", "PKM", "HK2", "MYC")
samples_plot <- c("rcc4_scrambled", "rcc4_mxi1_c1", "rcc4_mxi1_c2")
genes_plot <- c("KLF9", "CP")
# genes_plot <- c("KLF9", "HK2", "PFKP", "ENO2", "PKM")
# samples_plot <- c("rcc4_scrambled", "rcc4_klf9_c2", "rcc4_klf9_c3")
samples_plot <- c("rcc4_scrambled", "rcc4_klf9_c2")
sampletexts_plot <- c("sh-NC", "sh-KLF9")

# format expression data --------------------------------------------------
plot_data_long_df <- exp_df %>%
  filter(external_gene_name %in% genes_plot) %>%
  melt.data.table(measure.vars = colnames_value) %>%
  mutate(sample_text = gsub(x = variable, pattern = "sample\\.||_e1", replacement = "")) %>%
  filter(sample_text %in% samples_plot)
plot_data_long_df$sample_text2 <- mapvalues(x = plot_data_long_df$sample_text, from = samples_plot, to = sampletexts_plot)
plot_data_long_df$sample_text2 <- factor(x = plot_data_long_df$sample_text2, levels = sampletexts_plot)
plot_data_long_df$external_gene_name <- factor(x = plot_data_long_df$external_gene_name, levels = genes_plot)
p <- ggplot()
p <- p + geom_col(data = plot_data_long_df, mapping = aes(x = external_gene_name, y = value, fill = sample_text2), position=position_dodge(), color = "black")
p <- p + theme_classic()
p <- p + ylab(label = "% CPM to control")
p <- p + ggtitle(label = paste0("RNA-seq expression"))
# p <- p + theme(strip.background = element_rect(fill = NA),
#                panel.spacing = unit(0, "lines"))
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
p <- p + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank())
p
file2write <- paste0(dir_out, paste0(genes_plot, collapse = "_"), ".bulkRNA.CPM.", "pdf")
pdf(file2write, width = 3.5, height = 3, useDingbats = F)
print(p)
dev.off()
file2write <- paste0(dir_out, paste0(genes_plot, collapse = "_"), ".bulkRNA.CPM.", "png")
# png(file2write, width = 500, height = 450, res = 150)
png(file2write, width = 550, height = 450, res = 150)
print(p)
dev.off()


