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
## input degs
deg_all_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/unite_degs/unite_degs_for_C3L-00010_related_tumorcells_newcluster0_vs_others/20200818.v1/findallmarkers_C3L-00010_correlated_tumorcells_newcluster0vsothers.20200818.v1.tsv")

# make matrix for gene X cluster -------------------------------------------------------------
padj_mat <- dcast(data = deg_all_df, formula = gene ~ ident2, value.var = "p_val_adj")
avglogFC_mat <- dcast(data = deg_all_df, formula = gene ~ ident2, value.var = "avg_logFC")

# identify degs common to all the cluster comparisons and all significant ---------------------
## identify degs present in all the cluster comparisons
idx_allpresent <- (rowSums(is.na(padj_mat[,2:5])) == 0)
genes_allpresent <- padj_mat$gene[idx_allpresent]
genes_allpresent 
padj_nonna_mat <- padj_mat[idx_allpresent,]
## identify degs significant in all the cluster comparisons
idx_allsig <- (rowSums(padj_nonna_mat[2:5] > 0.05) == 0)
padj_allsig_mat <- padj_nonna_mat[idx_allsig,]
genes_allsig <- padj_allsig_mat$gene
genes_allsig

# identify degs with the same fold change directions ---------------------
avglogFC_allsig_mat <- avglogFC_mat[idx_allpresent,][idx_allsig,]
idx_samedirection <- (rowSums(avglogFC_allsig_mat[2:5] > 0) %in% c(0, 4))
genes_samedirection <- avglogFC_allsig_mat$gene[idx_samedirection]
genes_samedirection
avglogFC_samedirection_mat <- avglogFC_allsig_mat[idx_samedirection,]

# write output ------------------------------------------------------------
file2wirte <- paste0(dir_out, "shared_degs_C3L-00010_correlated_tumorcells_newcluster0vsothers", run_id, ".tsv")
write.table(x = avglogFC_samedirection_mat, file = file2wirte, quote = F, sep = "\t", row.names = F)
