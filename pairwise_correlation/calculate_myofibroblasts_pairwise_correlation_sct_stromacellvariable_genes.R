# Yige Wu @WashU March 2020
## running on local
## for calculating the correlation of average expression between tumor subclusters

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the variable gene list
variable_genes_df <- fread(input = "./Resources/Analysis_Results/findvariablefeatures/findvariablefeatures_cellgroup_stroma/20200811.v1/findvariablefeatures.cellgroup.Stroma.20200811.v1.tsv", data.table = F)
## input the average expression calculated (RNA)
avg.exp.mat <- fread(input = "./Resources/Analysis_Results/average_expression/averageexpression_sct_bycelltypeshorter_byaliquot_on_katmai/20200828.v1/averageexpression_SCT_bycelltypeshorter_byaliquot.31_aliquot_integration.20200828.v1.tsv", data.table = F)
avg.exp.mat <- avg.exp.mat %>%
  rename(gene = V1)

# format the column names to only aliquot id ------------------------------
## remove RNA from the column names
data_col_names <- colnames(avg.exp.mat)[-1]
data_col_names.changed <- str_split_fixed(string = data_col_names, pattern = "\\.", n = 2)[,2]
## rename the data frame
colnames(avg.exp.mat) <- c("gene", data_col_names.changed)
## remove the subcluster without NA as manual subcluster id
data_col_names.filtered <- data_col_names.changed[grepl(pattern = "_Myofibroblasts", x = data_col_names.changed)]
data_col_names.filtered
avg.exp.mat <- avg.exp.mat[, c("gene", data_col_names.filtered)]

# filter down to only variably expressed genes------------------------------------------------------------------
avg.exp.mat <- avg.exp.mat %>%
  filter(gene %in% variable_genes_df$gene)

# run and write pearson pairwise correlation ------------------------------------------------
pearson_coef <- cor(x = avg.exp.mat[, -1], use="complete.obs", method = "pearson")
pearson_coef
write.table(x = pearson_coef, file = paste0(dir_out, "pearson_coef.", "stromavaraible_genes.", "myofibroblasts.", "avg_exp.", "sct.", "byaliquot.", run_id, ".tsv"), quote = F, sep = "\t", row.names = T)

