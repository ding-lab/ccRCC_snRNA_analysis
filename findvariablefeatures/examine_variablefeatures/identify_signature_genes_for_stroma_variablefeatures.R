# Yige Wu @WashU Aug 2020
## reference: http://jokergoo.github.io/supplementary/ComplexHeatmap-supplementary1-4/supplS2_scRNASeq/supplS2_scRNAseq.html
## Signature genes are defined as a list of genes where each gene correlates to more than 20 genes with an absolute correlation larger than 0.5.

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
avgexp_df <- fread(input = "./Resources/Analysis_Results/average_expression/averageexpression_bycellgroup_byaliquot_on_katmai/20200804.v1/averageexpression_bycellgroup_byaliquot.30_aliquot_integration.20200804.v1.tsv", data.table = F)
avgexp_df <- avgexp_df %>%
  rename(gene = V1)

# format the column names to only aliquot id ------------------------------
## remove RNA from the column names
data_col_names <- colnames(avgexp_df)[-1]
data_col_names.changed <- str_split_fixed(string = data_col_names, pattern = "\\.", n = 2)[,2]
## rename the data frame
colnames(avgexp_df) <- c("gene", data_col_names.changed)
## remove the subcluster without NA as manual subcluster id
data_col_names.filtered <- data_col_names.changed[grepl(pattern = "_Stroma", x = data_col_names.changed)]
data_col_names.filtered
avgexp_df <- avgexp_df[, c("gene", data_col_names.filtered)]

# filter down to only variably expressed genes------------------------------------------------------------------
avgexp_df <- avgexp_df %>%
  filter(gene %in% variable_genes_df$gene)
avgexp_mat <- as.matrix(avgexp_df[,-1])
rownames(avgexp_mat) <- avgexp_df$gene

# identify signature genes ------------------------------------------------
## specify the cutoffs
cor_cutoff <- 0.5
n_cutoff <- 20
## calculate gene to gene correlation
dt = cor(t(avgexp_mat), method = "spearman")
dt %>% head()
diag(dt) = 0
dt[abs(dt) < cor_cutoff] = 0
dt[dt < 0] = -1
dt[dt > 0] = 1
i = colSums(abs(dt)) > n_cutoff
avgexp_signature_mat = avgexp_mat[i, ,drop = FALSE]
dim(avgexp_signature_mat)
