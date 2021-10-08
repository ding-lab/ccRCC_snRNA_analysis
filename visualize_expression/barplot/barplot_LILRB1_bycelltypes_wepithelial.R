# Yige Wu @WashU Sep 2021

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
avgexp_df <- fread(input = "./Resources/Analysis_Results/average_expression/avgexp_sct_data_bycelltypew_epithelial_katmai/20210907.v1/35_aliquot_merged.avgexp.SCT.data.Cell_group_w_epithelialcelltypes.20210907.v1.tsv", data.table = F)

# specify pairs to filter -------------------------------------------------
genes_filter <- c("LILRB1")

# format the column names to only aliquot id ------------------------------
## filtr the rows
plot_data_df <- avgexp_df %>%
  rename(gene = V1) %>%
  filter(gene %in% genes_filter) %>%
  melt() %>%
  mutate(cell_type = gsub(x = variable, pattern = "SCT\\.", replacement = "")) %>%
  filter(cell_type != "Unknown") %>%
  arrange(value)
plot_data_df$cell_type <- factor(x = plot_data_df$cell_type, levels = plot_data_df$cell_type)

# plot --------------------------------------------------------------------
p <- ggplot()
p <- p + geom_bar(data = plot_data_df, mapping = aes(x = value, y = cell_type), stat = "identity")
p <- p + xlab("Normalized  expression")
p <- p + ggtitle(paste0(genes_filter, " expression"))
p <- p + theme_classic(base_size = 14)
p <- p + theme(axis.text = element_text(color = "black", size = 14),
               axis.title.y = element_blank(),
               axis.title.x = element_text(size = 14),
               title = element_text(size = 14))
p

file2write <- paste0(dir_out, genes_filter, ".pdf")
pdf(file2write, width = 5, height = 4, useDingbats = F)
print(p)
dev.off()


