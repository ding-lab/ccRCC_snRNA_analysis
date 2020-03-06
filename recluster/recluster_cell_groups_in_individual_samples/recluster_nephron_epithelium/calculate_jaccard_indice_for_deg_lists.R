# Yige Wu @WashU March 2020
## fmake the matrix with the log2 fold change 

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the super table for degs
deg_sup_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/unite_degs/20200306.v1/FindAllMarkers.Wilcox.Pos.20200306.v1.tsv", data.table = F)
## crate the function to calculate jaccard index
jaccard <- function(M, user1, user2) {
  sums = rowSums(M[,c(user1, user2)])
  
  similarity = length(sums[sums==2])
  total = length(sums[sums==1]) + similarity
  
  similarity/total
}

# transform the data frame to matrix of binary values ---------------------
## create unique name for each tumor subcluster
deg_sup_df <- deg_sup_df %>%
  mutate(cluster_id = paste0(aliquot, "_", cluster)) %>%
  mutate(value = 1)
## transform the deg data frame to matrix of binary values
deg_binary_wide_df <- dcast(data = deg_sup_df, formula = gene ~ cluster_id, value.var = "value", fill = 0)
deg_binary_mat <- deg_binary_wide_df %>%
  select(-gene) %>%
  as.matrix()
rownames(deg_binary_mat) <- deg_binary_wide_df$gene
deg_binary_mat[1, 1:2]
## calculate jaccard indice from the matrix
jaccard(deg_binary_mat, "CPT0000870003_0", "CPT0000870003_1")

d <- sapply(colnames(deg_binary_mat), function(x, deg_binary_mat) sapply(names(deg_binary_mat), function(y, x, deg_binary_mat) jaccard(deg_binary_mat, y, x), deg_binary_mat = deg_binary_mat, x = x), deg_binary_mat = deg_binary_mat)
d
## create hclust object
hv <- hclust(as.dist(1-d))

## save hcluster object
