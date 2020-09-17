# Yige Wu @WashU Aug 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input not scaled average expression
avgexp_df <- fread(input = "./Resources/Analysis_Results/average_expression/averageexpression_sct_usescale_bycellgroup_on_katmai/20200904.v1/averageexpression_SCT_bycellgroup.detailed.31_aliquot_integration.20200904.v1.tsv", data.table = F)

# format the column names to only aliquot id ------------------------------
## filtr the rows
plot_data_df <- avgexp_df %>%
  rename(gene = V1)
## remove teh prefix from the column names
data_col_names <- colnames(plot_data_df)[-1]
data_col_names.changed <- str_split_fixed(string = data_col_names, pattern = "\\.", n = 2)[,2]
data_col_names.changed
## rename the data frame
colnames(plot_data_df) <- c("gene", data_col_names.changed)
plot_data_df <- plot_data_df[, c("gene", "Tumor.cells", "Immune", "Stroma", "Normal.epithelial.cells")]

# write table -------------------------------------------------------------
file2write <- paste0(dir_out, "formated.", "averageexpression.SCT.slotdata.bycellgroup.31_aliquot_integration.", run_id, ".tsv")
write.table(x = plot_data_df, file = file2write, quote = F, row.names = F, sep = "\t")

