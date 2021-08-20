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
library(lineup)

# input dependencies ------------------------------------------------------
## input the variable gene list from Lake et al
markers_lake_df <- fread(input = "./Resources/Analysis_Results/findvariablefeatures/findvariablefeatures_sct_Lake_bycluster_katmai/20210819.v1/findvariablefeatures.Lake_etal.20210819.v1.tsv", data.table = F)
## input the marker genes for our data
markers_ccrcc_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_by_celltype/unite_markers/unite_nephron_epithelial_celltype_markers/20210811.v1/ccRCC.Nephron_Epithelial_CellType_AlgorithmeticMarkers.20210811.v1.tsv")
## input the average expression for Lake et al

## input the average expression for our data
exp_df1 <- fread(input = "./Resources/Analysis_Results/average_expression/avgexp_sct_data_bycelltypew_epithelial_bysample_wPTreclustered_katmai/20210809.v1/35_samples_merged.avgexp.SCT.data.Cell_group_w_epithelialcelltypes20210809.v1.tsv", data.table = F)
exp_df2 <- fread(data.table = F, input = "./Resources/Analysis_Results/average_expression/avgexp_Lake_sct_data_bycluster_katmai/20210819.v1/AverageExpression.SCT.dataLake_etal.20210819.v1.tsv")
# process the genes to use to do pairwise correlation ---------------------
markers_process_df <- markers_ccrcc_df %>%
  filter(Gene %in% markers_lake_df$gene)
table(markers_process_df$cell_type)
# Distal convoluted tubule       Intercalated cells            Loop of Henle                Podocytes          Principle cells          Proximal tubule 
# 47                       79                       63                       61                       48                      123 
# Tumor cells 
# 111 

# format expression data --------------------------------------------------
exp_mat1 <- as.matrix(exp_df1[,-1]); rownames(exp_mat1) <- exp_df1$V1
exp_filtered_mat1 <- exp_mat1[markers_process_df$Gene,]
exp_mat2 <- as.matrix(exp_df2[,-1]); rownames(exp_mat2) <- exp_df2$V1
exp_filtered_mat2 <- exp_mat2[markers_process_df$Gene,]


# calculate correalation --------------------------------------------------
cor_result <- cor(x = exp_filtered_mat1, y = exp_filtered_mat2, method = "spearman")
file2write <- paste0(dir_out, "ccRCC_vs_Lake_etal.", "spearman.", "sct.data.", "celltypemarkers_overlap_variablegenes.", "tsv")
write.table(x = cor_result, file = file2write, quote = F, sep = "\t", row.names = T)
cor_result <- cor(x = exp_filtered_mat1, y = exp_filtered_mat2, method = "pearson")
file2write <- paste0(dir_out, "ccRCC_vs_Lake_etal.", "pearson", "sct.data.", "celltypemarkers_overlap_variablegenes.", "tsv")
write.table(x = cor_result, file = file2write, quote = F, sep = "\t", row.names = T)


