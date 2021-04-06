# Yige Wu @WashU Mar 2021

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
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/pbrm1_bap1_vs_non_mutants/summarize_degs/summarize_PBRM1_BAP1_DEGs/20210405.v1/BAP1_PBRM1_DEGs.Num_samples.20210405.v1.tsv")

# make plot data ----------------------------------------------------------
plot_data_df <- deg_df %>%
  select(genesymbol_deg, PBRM1_vs_BAP1_Mutants_snRNA, PBRM1_vs_NonMutants_snRNA, PBRM1_Tumorcells_vs_PTcells_snRNA) %>%
  melt(id.var = "genesymbol_deg") %>%
  filter(value != "Inconsistent") %>%
  group_by(variable, value) %>%
  summarize(number_degs = n()) %>%
  mutate(Direction = value)
## make colors
color_red <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
color_blue <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[2]

# plot --------------------------------------------------------------------
p <- ggplot()
p <- p + geom_bar(data = plot_data_df, mapping = aes(x = variable, y = number_degs, fill = Direction), stat = "identity")
# p <- p + theme(axis.text.x = element_text(angle = 15))
p <- p + coord_flip()
p <- p + scale_fill_manual(values = c("Up" = color_red, "Down" = color_blue))
p <- p + theme_classic()
p <- p + theme(axis.text.y = element_text(size = 10), axis.title.y = element_blank())
p
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "PBRM1.", "Number_DEGs_byDirection.", "png")
png(file2write, width = 850, height = 200, res = 150)
print(p)
dev.off()


