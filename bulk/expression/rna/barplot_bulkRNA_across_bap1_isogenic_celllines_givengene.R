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
  mutate(sample_text = gsub(pattern = "X", replacement = "", x = variable)) %>%
  mutate(log2TPM = log2(value+1))

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
file2write <- paste0(dir_out, gene_plot, ".bulkRNA.", "pdf")
pdf(file2write, width = 3, height = 3.5, useDingbats = F)
print(p)
dev.off()

