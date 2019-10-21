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
int_run_id <- "20191015.v1"
int_seurat_obj <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/analysis_results/integration/integrate_seurat_objects/20191015.v1/Renal_Integrated.20191015.v1.RDS")
aliquot_ids <- unique(int_seurat_obj@meta.data$orig.ident)

# get cluster-average of sn gene expression ----------------------------------
## in the integrated assay, there are 2000 genes because the integrated used the variable genes (maybe as anchor) and it only kept the top 2000 most variably expressed genes
## however, the PDCD1/PD1 is not within the 2000 genes

## for how averageexpression is calculated, refer to
## https://github.com/satijalab/seurat/issues/705
## https://github.com/satijalab/seurat/blob/aafba71530ae5bfd602eea380b8b45f0062580ad/R/utilities.R#L227-L277
sn_ave_scaled_exp_tab <- AverageExpression(object = int_seurat_obj, assays = "integrated", use.scale = T, add.ident = "orig.ident")
sn_ave_scaled_exp_tab$integrated %>% nrow()
sn_ave_scaled_exp_tab <- sn_ave_scaled_exp_tab$integrated
sn_ave_scaled_exp_tab %>% head()
write.table(x = sn_ave_scaled_exp_tab, file = paste0(dir_out, "Average_Scaled_Data_for_Integrated.", int_run_id,".Run.", run_id, ".tsv"), quote = F, row.names = T, col.names = T, sep = "\t")

sn_ave_counts_tab <- AverageExpression(object = int_seurat_obj, assays = "RNA", slot = "counts", add.ident = "orig.ident")
sn_ave_counts_tab$RNA %>% nrow()
sn_ave_counts_tab <- sn_ave_counts_tab$RNA
sn_ave_counts_tab %>% head()
write.table(x = sn_ave_counts_tab, file = paste0(dir_out, "Average_RNA_Counts_for_Integrated.", int_run_id, ".Run.", run_id,".tsv"), quote = F, row.names = T, col.names = T, sep = "\t")

# input protein data ------------------------------------------------------
pro_tab <- loadParseProteomicsData(expression_type  = "PRO", sample_type = "tumor")
pro_mat <- pro_tab
rownames(pro_mat) <- pro_tab$Gene
pro_mat$Gene <- NULL
pro_mat <- as.matrix(pro_mat)

# reformat protein data colnames to aliquot ids ---------------------------
## input meta data
meta_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/analysis_results/sample_info/make_meta_data/meta_data.20190924.v1.tsv", data.table = F)

meta_tab_in_use <- meta_tab %>%
  dplyr::filter(meta_tab$Specimen.ID.snRNA %in% aliquot_ids)
rownames(meta_tab_in_use) <- meta_tab_in_use$Case.ID

pro_mat <- pro_mat[,meta_tab_in_use$Case.ID]
pro_mat
colnames(pro_mat) <- meta_tab_in_use$Specimen.ID.snRNA

# run spearman correlation with the snRNA scaled data per gene---------------------------
## define the genes to test
genes2test <- intersect(rownames(sn_ave_scaled_exp_tab), pro_tab$Gene)
length(genes2test)
## 1161 genes

"CD274" %in% genes2test
"CD38" %in% genes2test
"POSTN" %in% genes2test
"PDCD1LG2" %in% genes2test

## run correlation by cluster by gene
spearman.test.result.df.sup <- NULL
for (cluster_tmp in unique(int_seurat_obj@meta.data$seurat_clusters)) {
  sn_exp_cluster_tmp <- sn_ave_scaled_exp_tab[, paste0(cluster_tmp, "_", meta_tab_in_use$Specimen.ID.snRNA)]
  
  spearman.test.result.df <- ldply(genes2test, function(g, exp1, exp2) {
    exp1_tmp <- exp1[g,] %>% unlist()
    exp2_tmp <- exp2[g,] %>% unlist()
    if (sum(complete.cases(data.frame(a = exp1_tmp, b = exp2_tmp))) < 4) {
      spearman.test.result.tmp <- c(gene_symbol = g, p.value = NA, rho = NA)
    } else {
      spearman.test.tmp <- cor.test(x = exp1_tmp, y = exp2_tmp, method = "spearman")
      spearman.test.result.tmp <- c(gene_symbol = g, p.value = spearman.test.tmp$p.value, spearman.test.tmp$estimate)
    }
    return(spearman.test.result.tmp)
  }, exp1 = sn_exp_cluster_tmp, exp2 = pro_mat)
  spearman.test.result.df$cluster <- cluster_tmp
  spearman.test.result.df.sup <- rbind(spearman.test.result.df.sup, spearman.test.result.df)
}
spearman.test.result.df.sup <- spearman.test.result.df.sup %>%
  arrange(p.value)
write.table(x = spearman.test.result.df.sup, 
            file = paste0(dir_out, "Spearman_Correlation_btw_Average_Scaled_snRNA_data_w_bulk_Protein_for_Integrated.", int_run_id,".Run.", run_id, ".tsv"), 
            quote = F, row.names = F, col.names = T, sep = "\t")


# run spearman correlation with the snRNA count data per gene---------------------------
## define the genes to test
genes2test <- intersect(rownames(sn_ave_counts_tab), pro_tab$Gene)
length(genes2test)
## 1161 genes

"CD274" %in% genes2test
"CD38" %in% genes2test
"POSTN" %in% genes2test
"PDCD1LG2" %in% genes2test

## run correlation by cluster by gene
spearman.test.result.df.sup <- NULL
for (cluster_tmp in unique(int_seurat_obj@meta.data$seurat_clusters)) {
  sn_exp_cluster_tmp <- sn_ave_counts_tab[, paste0(cluster_tmp, "_", meta_tab_in_use$Specimen.ID.snRNA)]
  
  spearman.test.result.df <- ldply(genes2test, function(g, exp1, exp2) {
    exp1_tmp <- exp1[g,] %>% unlist()
    exp2_tmp <- exp2[g,] %>% unlist()
    if (sum(complete.cases(data.frame(a = exp1_tmp, b = exp2_tmp))) < 4) {
      spearman.test.result.tmp <- c(gene_symbol = g, p.value = NA, rho = NA)
    } else {
      spearman.test.tmp <- cor.test(x = exp1_tmp, y = exp2_tmp, method = "spearman")
      spearman.test.result.tmp <- c(gene_symbol = g, p.value = spearman.test.tmp$p.value, spearman.test.tmp$estimate)
    }
    return(spearman.test.result.tmp)
  }, exp1 = sn_exp_cluster_tmp, exp2 = pro_mat)
  spearman.test.result.df$cluster <- cluster_tmp
  spearman.test.result.df.sup <- rbind(spearman.test.result.df.sup, spearman.test.result.df)
}
spearman.test.result.df.sup <- spearman.test.result.df.sup %>%
  arrange(p.value)
write.table(x = spearman.test.result.df.sup, 
            file = paste0(dir_out, "Spearman_Correlation_btw_Average_Count_snRNA_data_w_bulk_Protein_for_Integrated.", int_run_id,".Run.", run_id, ".tsv"), 
            quote = F, row.names = F, col.names = T, sep = "\t")
