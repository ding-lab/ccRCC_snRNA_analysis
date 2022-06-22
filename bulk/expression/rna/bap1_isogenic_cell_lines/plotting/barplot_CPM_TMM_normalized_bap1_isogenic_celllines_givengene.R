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
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
exp_df <- fread(input = "./Resources/Analysis_Results/bulk/expression/rna/other/normalization/get_edgeR_TMM_normalized_counts_all_cell_lines/20220606.v1/CPM.TMM_normalized.All_Cell_Lines.20220606.v1.tsv", data.table = F)

# specify parameters ------------------------------------------------------
gene_plot <- "BAP1"
gene_plot <- "CES3"
# gene_plot <- "PTPRJ"

# format expression data --------------------------------------------------
colnames_id <- colnames(exp_df)[!grepl(x = colnames(exp_df), pattern = "sample")]
plot_data_long_df <- exp_df %>%
  filter(external_gene_name %in% gene_plot) %>%
  melt(id.vars = colnames_id) %>%
  # filter(variable %in% c("sample.786_o.sglacz_e1", "sample.786_o.sgbap1_e1", "sample.skrc42.bap1_e1", "sample.skrc42.emptyvector_e1")) %>%
  filter(variable %in% c("sample.skrc42.bap1_e1", "sample.skrc42.emptyvector_e1")) %>%
  mutate(parental_line = str_split_fixed(string = variable, pattern = "\\.", n = 3)[,2]) %>%
  mutate(parental_line_text = ifelse(parental_line == "786_o", "786-O", "SKRC-42")) %>%
  mutate(BAP1_status = ifelse(variable %in% c("sample.786_o.sgbap1_e1", "sample.skrc42.emptyvector_e1"), "BAP1 null", "BAP1 wt"))
plot_data_long_df$sample_text <- as.vector(plot_data_long_df$variable)
plot_data_long_df$sample_text[plot_data_long_df$sample_text == "sample.786_o.sgbap1_e1"] <- "sgBAP1"
plot_data_long_df$sample_text[plot_data_long_df$sample_text == "sample.786_o.sglacz_e1"] <- "sgLacZ"
plot_data_long_df$sample_text[plot_data_long_df$sample_text == "sample.skrc42.bap1_e1"] <- "BAP1"
plot_data_long_df$sample_text[plot_data_long_df$sample_text == "sample.skrc42.emptyvector_e1"] <- "control"
plot_data_long_df$sample_text <- factor(x = plot_data_long_df$sample_text, levels = c("sgLacZ", "sgBAP1", "BAP1", "control"))

# make barplot ------------------------------------------------------------

p <- ggplot()
p <- p + geom_col(data = plot_data_long_df, mapping = aes(x = sample_text, y = value, fill = BAP1_status), color = "black")
# p <- p + facet_wrap(.~parental_line_text, scales = "free_x")
p <- p + theme_classic()
p <- p + ylab(label = paste0(gene_plot, " expression (TPM)"))
p <- p + theme(strip.background = element_rect(fill = NA, color = NA),
               panel.spacing = unit(0, "lines"))
p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, color = "black"), axis.text.y = element_text(color = "black"))
p <- p + theme(axis.title.x = element_blank(), axis.ticks.x = element_blank())
p
file2write <- paste0(dir_out, gene_plot, ".bulkRNA.TPM.", "pdf")
# pdf(file2write, width = 3, height = 2.5, useDingbats = F)
pdf(file2write, width = 2.1, height = 2, useDingbats = F)
print(p)
dev.off()
file2write <- paste0(dir_out, gene_plot, ".bulkRNA.TPM.", "png")
# png(file2write, width = 425, height = 300, res = 150)
png(file2write, width = 325, height = 300, res = 150)
print(p)
dev.off()

