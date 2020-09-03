# Yige Wu @WashU Aug 2020

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
genes_celltypemarker_df <- fread(input = "./Resources/Knowledge/Kidney_Markers/Gene2CellType_Tab.20200831.v1.tsv", data.table = F)
## input the average expression calculated (RNA)
avg.exp.mat <- fread(input = "./Resources/Analysis_Results/average_expression/averageexpression_sct_bycelltypeshorter_byaliquot_on_katmai/20200828.v1/averageexpression_SCT_bycelltypeshorter_byaliquot.31_aliquot_integration.20200828.v1.tsv", data.table = F)
## input meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")

# specify the genes to show -----------------------------------------------
genes_plot <- genes_celltypemarker_df$Gene[genes_celltypemarker_df$Cell_Type1 %in% c("Fibroblasts", "Pericyte", "Vascular smooth muscle cells and pericytes")]

# format the column names to only aliquot id ------------------------------
plot_data_df <- avg.exp.mat %>%
  rename(gene = V1) %>%
  filter(gene %in% genes_plot)
## remove RNA from the column names
data_col_names <- colnames(plot_data_df)[-1]
## get only myofibroblasts
data_col_names.filtered <- data_col_names[grepl(pattern = "_Myofibroblasts", x = data_col_names)]
data_col_names.filtered
## filter data matrix
plot_data_mat <- as.matrix(plot_data_df[, c(data_col_names.filtered)])
## change column names
data_col_names.changed <- str_split_fixed(string = data_col_names.filtered, pattern = "\\.", n = 2)[,2]
## rename the data frame
colnames(plot_data_mat) <- data_col_names.changed
rownames(plot_data_mat) <- plot_data_df$gene
plot_data_mat %>% head()

# write outpu -------------------------------------------------------------
file2write <- paste0(dir_out, "myofibroblasts_markergenes.", "myofibroblasts", "bysample.", ".tsv")
write.table(x = plot_data_mat, file = file2write, quote = F, sep = "\t", row.names = T)
