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
avgexp_df <- fread(input = "./Resources/Analysis_Results/average_expression/averageexpression_sct_usescale_bycellgroup3_byaliquot_on_katmai/20200917.v1/avgexp.SCT.bycellgroup3.byaliquot.31_aliquot_integration.20200917.v1.tsv", data.table = F)
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")
# specify the gene to plot ------------------------------------------------
gene_plot <- "VIM"

# format expression data --------------------------------------------------
plot_data_df <- avgexp_df %>%
  filter(V1 == gene_plot) %>%
  melt() %>%
  mutate(id_aliquot_cellgroup3 = gsub(x = variable, pattern = "SCT.", replacement = "")) %>%
  mutate(aliquot = str_split_fixed(string = id_aliquot_cellgroup3, pattern = "_", n = 2)[,1]) %>%
  mutate(cellgroup3 = str_split_fixed(string = id_aliquot_cellgroup3, pattern = "_", n = 2)[,2])
plot_data_df <- plot_data_df %>%
  filter(cellgroup3 == "Nephron_Epithelium")
plot_data_df$id_aliquot_wu <- mapvalues(x = plot_data_df$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))


# write output ------------------------------------------------------------
p

