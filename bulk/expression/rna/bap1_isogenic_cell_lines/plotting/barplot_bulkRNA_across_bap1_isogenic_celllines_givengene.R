# Yige Wu @WashU Dec 2020

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
## input the average expression calculated (SCT)
exp_df <- fread(input = "./Resources/Analysis_Results/bulk/expression/rna/unite_bap1_isogenic_cellline_mrna_kallisto/20211221.v1/BAP1_isogenic_celllines.Kallisto.gene_level.TPM.20211221.v1.tsv", data.table = F)

# specify parameters ------------------------------------------------------
gene_plot <- "CES3"

# format expression data --------------------------------------------------
plot_data_long_df <- exp_df %>%
  filter(gene_name %in% gene_plot) %>%
  melt()
plot_data_long_df <- plot_data_long_df %>%
  mutate(log2TPM = log2(value+1))
plot_data_long_df$sample_text <- as.vector(plot_data_long_df$variable)
plot_data_long_df$sample_text[plot_data_long_df$sample_text == "X786.o_sgbap1"] <- "786O_BAP1mt"
plot_data_long_df$sample_text[plot_data_long_df$sample_text == "X786.o_sglacz"] <- "786O_BAP1wt"
plot_data_long_df$sample_text[plot_data_long_df$sample_text == "skrc42_bap1"] <- "SKRC42_BAPwt"
plot_data_long_df$sample_text[plot_data_long_df$sample_text == "skrc42_emptyvector"] <- "SKRC42_BAPmt"
plot_data_long_df$sample_text <- factor(x = plot_data_long_df$sample_text, levels = c("SKRC42_BAPwt", "SKRC42_BAPmt", "786O_BAP1wt", "786O_BAP1mt"))

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

