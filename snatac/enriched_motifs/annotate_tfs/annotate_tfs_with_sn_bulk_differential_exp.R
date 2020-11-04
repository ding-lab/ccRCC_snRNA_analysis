# Yige Wu @WashU Oct 2020

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
## input the DEG-TF matrix
deg2motif_wide_df <- fread(data.table = F, input = "./Resources/snRNA_Processed_Data/Differentially_Expressed_Genes/DEGs_with_TFs_inDARs.tsv")
## input the deg table
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/findallmarker_wilcox_tumor_vs_pt_on_katmai/20200929.v1/findallmarkers_wilcox_tumorcells_vs_pt.20200929.v1.tsv")
## input the differentially expression of bulk 
bulk_rna_de_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/Tumor_normal_bulk_RNA.txt")
bulk_pro_de_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/Tumor_normal_bulk_Protein.txt")

# get unique motifs and TF names ------------------------------------------
## get motif names
motifnames_vec <- colnames(deg2motif_wide_df)[-1]
## format the motif names to TF names
### duplicate rows with multiple gene symbol
idx_rep <- sapply(X = motifnames_vec, FUN = function(x) {
  vec_genes <- str_split(string = x, pattern = "\\:")[[1]]
  vec_genes <- gsub(x = vec_genes, pattern = "(var.2)", replacement = "")
  vec_genes <- vec_genes[vec_genes != ""]
  length_genes <- length(vec_genes)
  return(length_genes)
})
motif2genesymbol_df <- data.frame(motifname = rep(motifnames_vec, idx_rep))
### get the gene symbols of the TFs
genesymbols_tf <- sapply(X = motifnames_vec, FUN = function(x) {
  vec_genes <- str_split(string = x, pattern = "\\:")[[1]]
  vec_genes <- gsub(x = vec_genes, pattern = "\\(var.2\\)", replacement = "")
  vec_genes <- vec_genes[vec_genes != ""]
  return(vec_genes)
})
motif2genesymbol_df$genesymbol_tf <- unlist(genesymbols_tf)

# merge with snRNA DEGs ---------------------------------------------------
tfgene2de_df <- data.frame(genesymbol_tf = unique(motif2genesymbol_df$genesymbol_tf))
## merge with snRNA
colnames(deg_df) <- paste0(colnames(deg_df), ".snrna")
tfgene2de_df <- merge(x = tfgene2de_df, y = deg_df, 
                      by.x = c("genesymbol_tf"), by.y = c("row_name.snrna"), all.x = T)
## merge with bulk mRNA
colnames(bulk_rna_de_df) <- paste0(colnames(bulk_rna_de_df), ".bulk.rna")
tfgene2de_df <- merge(x = tfgene2de_df, y = bulk_rna_de_df, 
                      by.x = c("genesymbol_tf"), by.y = c("Gene.bulk.rna"), all.x = T)
## merge with bulk protein
colnames(bulk_pro_de_df) <- paste0(colnames(bulk_pro_de_df), ".bulk.pro")
tfgene2de_df <- merge(x = tfgene2de_df, y = bulk_pro_de_df, 
                      by.x = c("genesymbol_tf"), by.y = c("Gene.bulk.pro"), all.x = T)
## sort
tfgene2de_df <- tfgene2de_df %>%
  arrange(FDR.bulk.pro)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Top_TF_DE_snRNA_and_bulk", ".tsv")
write.table(x = tfgene2de_df, file = file2write, quote = F, sep = "\t", row.names = F)



