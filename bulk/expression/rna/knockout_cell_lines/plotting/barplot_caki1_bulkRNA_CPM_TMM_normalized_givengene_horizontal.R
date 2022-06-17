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
exp_df <- fread(input = "./Resources/Analysis_Results/bulk/expression/rna/knockout_cell_lines/normalization/get_edgeR_TMM_normalized_counts_selected_cell_lines/20220524.v1/CPM.TMM_normalized.Knockout_Cell_Lines.20220524.v1.tsv", data.table = F)
# exp_df <- fread(input = "./Resources/Analysis_Results/bulk/expression/rna/knockout_cell_lines/divide_knockout_cellline_TPM_TMMnormalized_byscrambled/20220516.v1/Knockout_celllines.Kallisto.gene_level.CPM.TMMnormalized.DividedByControl.20220516.v1.tsv", data.table = F)

# specify parameters ------------------------------------------------------
colnames_value <- colnames(exp_df)[grepl(pattern = "sample", x = colnames(exp_df))]
exp_df <- as.data.table(exp_df)
# genes_plot <- c("KLF9", "CP")
# genes_plot <- c("CP")
genes_plot <- c("COL4A1")
samples_plot <- c("caki_1_control_e1", "dr_caki_1_rna", "caki_1_cp_c2_e1", "caki_1_cp_c1_e1")
sampletext_plot <- c("caki1_nt1", "caki1_nt2", "caki1_cp_c2", "caki1_cp_c1")
sampletext_plot <- c("sh-NT1", "sh-NT2", "sh-CP-C2", "sh-CP-C1")

# make colors -------------------------------------------------------------
color_nt <- RColorBrewer::brewer.pal(n = 6, name = "Set3")[6]
color_cp <- RColorBrewer::brewer.pal(n = 6, name = "Set3")[5]
colors_bysample <- c(color_nt, color_nt, color_cp, color_cp)
names(colors_bysample) <- c("sh-NT1", "sh-NT2", "sh-CP-C2", "sh-CP-C1")

# format expression data --------------------------------------------------
plot_data_long_df <- exp_df %>%
  filter(external_gene_name %in% genes_plot) %>%
  melt.data.table(measure.vars = colnames_value) %>%
  mutate(sample = gsub(x = variable, pattern = "sample\\.", replacement = "")) %>%
  filter(sample %in% samples_plot)
plot_data_long_df$sample_text <- mapvalues(x = plot_data_long_df$sample, from = samples_plot, to = sampletext_plot)
plot_data_long_df$sample_text <- factor(x = plot_data_long_df$sample_text, levels = rev(sampletext_plot))

p <- ggplot()
p <- p + geom_col(data = plot_data_long_df, mapping = aes(x = value, y = sample_text, fill = sample_text), position=position_dodge())
p <- p + scale_fill_manual(values = colors_bysample)
p <- p + theme_classic()
p <- p + xlab(label = paste0(genes_plot, " gene CPM"))
# p <- p + ggtitle(label = paste0(paste0(genes_plot, collapse = " & "), " expression")) # , subtitle = "by RNA-seq"
# p <- p + theme(strip.background = element_rect(fill = NA),
#                panel.spacing = unit(0, "lines"))
p <- p + theme(axis.text = element_text(size = 15, color = "black"))
p <- p + theme(axis.title.y = element_blank(), axis.title.x = element_text(size = 15),
               axis.ticks.y = element_blank(), legend.position = "none")
p
file2write <- paste0(dir_out, paste0(genes_plot, collapse = "_"), ".bulkRNA.CPM.", "pdf")
pdf(file2write, width = 2.5, height = 1.25, useDingbats = F)
print(p)
dev.off()
# file2write <- paste0(dir_out, paste0(genes_plot, collapse = "_"), ".bulkRNA.CPM.", "png")
# png(file2write, width = 300, height = 450, res = 150)
# print(p)
# dev.off()

## when genes_plot == "CP
exp_shRNA <- sum(plot_data_long_df$value[plot_data_long_df$sample %in% c("caki_1_cp_c2_e1", "caki_1_cp_c1_e1")])
exp_NT <- sum(plot_data_long_df$value[plot_data_long_df$sample %in% c("caki_1_control_e1", "dr_caki_1_rna")])
exp_shRNA/exp_NT

