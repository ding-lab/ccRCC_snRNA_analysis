# Yige Wu @WashU Oct 2019
## for testing spearman correlation between cluster-average sn gene expression with bulk protein/phosphorylation level

# source ------------------------------------------------------------------
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")

# set run id  ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)


# set output directory ----------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# input integrated and cell type assigned seurat object ------------------------------------------
int_seurat_obj <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/analysis_results/integration/integrate_seurat_objects/20191002.v1/renal_integrated20191002.v1.RDS")

# get cluster-average of sn gene expression ----------------------------------
## in the integrated assay, there are 2000 genes because the integrated used the variable genes (maybe as anchor) and it only kept the top 2000 most variably expressed genes
## however, the PDCD1/PD1 is not within the 2000 genes

## for how averageexpression is calculated, refer to
## https://github.com/satijalab/seurat/issues/705
## https://github.com/satijalab/seurat/blob/aafba71530ae5bfd602eea380b8b45f0062580ad/R/utilities.R#L227-L277
sn_ave_scaled_exp_tab <- AverageExpression(object = int_seurat_obj, assays = "integrated", use.scale = T)
sn_ave_scaled_exp_tab$integrated %>% nrow()
sn_ave_scaled_exp_tab <- sn_ave_scaled_exp_tab$integrated
write.table(x = sn_ave_scaled_exp_tab, file = paste0(dir_out, ".Average_Scaled_Data_for_Integrated_3FACS.", "20191002.v1",".Run.", run_id, ".tsv"), quote = F, row.names = T, col.names = T, sep = "\t")

sn_ave_counts_tab <- AverageExpression(object = int_seurat_obj, assays = "RNA", slot = "counts")
sn_ave_counts_tab$RNA %>% nrow()
sn_ave_counts_tab <- sn_ave_counts_tab$RNA
sn_ave_counts_tab %>% head()

write.table(x = sn_ave_counts_tab, file = paste0(dir_out, ".Average_RNA_Counts_for_Integrated_3FACS.", "20191002.v1", ".Run.", run_id,".tsv"), quote = F, row.names = T, col.names = T, sep = "\t")

# input protein data ------------------------------------------------------
pro_tab <- loadParseProteomicsData(expression_type  = "PRO", sample_type = "tumor")
pro_mat <- pro_tab
rownames(pro_mat) <- pro_tab$Gene
pro_mat$Gene <- NULL
pro_mat <- as.matrix(pro_mat)
# pho_tab <- loadParseProteomicsData(expression_type  = "PHO", sample_type = "tumor")


# run spearman correlation with the scaled data ---------------------------
## define the genes to test
genes2test <- intersect(rownames(sn_ave_scaled_exp_tab), pro_tab$Gene)
length(genes2test)

## 1161 genes
"CD274" %in% genes2test
"CD38" %in% genes2test
"POSTN" %in% genes2test
"PDCD1LG2" %in% genes2test

## 
sn_ave_scaled_exp_tab[genes2test[1],]
lapply(head(genes2test), function(g, exp1, exp2) {
  exp1_tmp <- exp1[g,]
  exp2_tmp <- exp2[g,]
  
}, exp1 = sn_ave_scaled_exp_tab, exp2 = pro_mat)


# calculate median protein level ------------------------------------------
pro_tab <- pro_tab[pro_tab$Gene %in% genes2test,]

# pho_tab <- pho_tab[pho_tab$Gene %in% genes2test,]

