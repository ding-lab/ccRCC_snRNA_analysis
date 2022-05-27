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
exp_df <- fread(input = "./Resources/Analysis_Results/bulk/expression/rna/unite_knockout_cellline_mrna_kallisto/20220510.v1/Knockout_celllines.Kallisto.gene_level.TPM.20220510.v1.tsv", data.table = F)

# specify parameters ------------------------------------------------------
gene_plot <- "MXI1"
for (gene_plot in c("MXI1", "KLF9", "CP", "PCSK6", "HK2", "PKM", "PFKP", "ENO2", "MYC")) {
  # format expression data --------------------------------------------------
  plot_data_long_df <- exp_df %>%
    filter(gene_name %in% gene_plot) %>%
    melt()
  plot_data_long_df <- plot_data_long_df %>%
    mutate(log2TPM = log2(value+1))
  plot_data_long_df$sample_text <- as.vector(plot_data_long_df$variable)
  # plot_data_long_df$sample_text <- factor(x = plot_data_long_df$sample_text, levels = c("SKRC42_BAPwt", "SKRC42_BAPmt", "786O_BAP1wt", "786O_BAP1mt"))
  
  # make barplot ------------------------------------------------------------
  p <- ggplot()
  p <- p + geom_col(data = plot_data_long_df, mapping = aes(x = sample_text, y = log2TPM))
  p <- p + theme_classic()
  p <- p + ylab(label = "log2(TPM+1)")
  p <- p + ggtitle(label = paste0(gene_plot, " RNA expression"))
  # p <- p + theme(strip.background = element_rect(fill = NA),
  #                panel.spacing = unit(0, "lines"))
  p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  p <- p + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank())
  p
  file2write <- paste0(dir_out, gene_plot, ".bulkRNA.log2TPM.", "pdf")
  pdf(file2write, width = 2.5, height = 3, useDingbats = F)
  print(p)
  dev.off()
  file2write <- paste0(dir_out, gene_plot, ".bulkRNA.log2TPM.", "png")
  png(file2write, width = 400, height = 450, res = 150)
  print(p)
  dev.off()
  
  p <- ggplot()
  p <- p + geom_col(data = plot_data_long_df, mapping = aes(x = sample_text, y = value))
  p <- p + theme_classic()
  p <- p + ylab(label = "TPM")
  p <- p + ggtitle(label = paste0(gene_plot, " RNA expression"))
  # p <- p + theme(strip.background = element_rect(fill = NA),
  #                panel.spacing = unit(0, "lines"))
  p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
  p <- p + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank())
  p
  file2write <- paste0(dir_out, gene_plot, ".bulkRNA.TPM.", "pdf")
  pdf(file2write, width = 2.5, height = 3, useDingbats = F)
  print(p)
  dev.off()
  file2write <- paste0(dir_out, gene_plot, ".bulkRNA.TPM.", "png")
  png(file2write, width = 400, height = 450, res = 150)
  print(p)
  dev.off()
}

