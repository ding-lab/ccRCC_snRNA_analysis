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
## input methylation
methyl_df <- fread(data.table = F, input = "~/Box/CPTAC_ccRCC/Data_Freeze_1.1/CPTAC_ccRCC_Discovery/DNA_methylation/CPTAC_ccRCC_discovery_tumor_methylation_betavalue_probe_level_v1.0.tsv")
## input methylation annotation
methyl_anno_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Methylation/EPIC.hg38.manifest.tsv")
## input gene expression
rna_df <- fread(data.table = F, input = "~/Box/CPTAC_ccRCC/Data_Freeze_1.1/CPTAC_ccRCC_Discovery/Gene_expression/CPTAC_ccRCC_discovery_tumor_mRNA_FPKM_UQ_log2_v1.0.tsv")
## input clear-cell RCC classification
histology_df <- fread(data.table = F, input = "~/Box/CPTAC_ccRCC/Data_Freeze_1.1/CPTAC_ccRCC_Discovery/Case_ID/CPTAC_ccRCC_discovery_caseID_v1.0.tsv")

# filter methylation ------------------------------------------------------
nrow(methyl_df) ## 917437 rows
## some of probes are duplicated, and the values are the same across samples
methyl_df$Locus[!(methyl_df$Locus %in% methyl_anno_df$probeID)]
methyl_df$Locus[duplicated(methyl_df$Locus)]
methyl_df[methyl_df$Locus == "cg00000714",]
## only use probes with gene annotated
probeids_w_gene <- methyl_anno_df$probeID[!is.na(methyl_anno_df$gene)]
length(probeids_w_gene) ## 706421
## filter
methyl_df <- methyl_df[methyl_df$Locus %in% probeids_w_gene,]
# methyl_df <- unique(methyl_df)
methyl_df <- methyl_df[!duplicated(methyl_df$Locus),]

# make probe2gene map -----------------------------------------------------
probe2genes_df <- methyl_anno_df %>%
  filter(probeID %in% methyl_df$Locus) %>%
  select(probeID, gene_HGNC)
idx_rep <- sapply(probe2genes_df$gene_HGNC, function(x) {
  genes_vec <- str_split(string = x, pattern = ";")[[1]]
  genes_len <- length(genes_vec)
  return(genes_len)
})
probe2gene_df <- probe2genes_df[rep(1:nrow(probe2genes_df), idx_rep),]
probe2gene_df$gene_HGNC <- unlist(sapply(probe2genes_df$gene_HGNC, function(x) {
  genes_vec <- str_split(string = x, pattern = ";")[[1]]
  return(genes_vec)
}))
rm(methyl_anno_df)
## write output
file2write <- paste0(dir_out, "Probe2Genes_HGNC.", run_id, ".tsv")
write.table(x = probe2genes_df, file = file2write, quote = F, sep = "\t", row.names = F)
rm(probe2genes_df)
file2write <- paste0(dir_out, "Probe2Gene_HGNC.", run_id, ".tsv")
write.table(x = probe2gene_df, file = file2write, quote = F, sep = "\t", row.names = F)

# process the methylation and RNA matrix again ---------------------------------------------------------
## figure out how many genes will be tested
genes2test <- rna_df$gene_name[rna_df$gene_name %in% probe2gene_df$gene_HGNC]
length(genes2test) ## 46206
## figure out which cases will be tested
caseids2test <- histology_df$CASE_ID[histology_df$Histologic_Type == "Clear cell renal cell carcinoma"]
length(caseids2test) ## 103
## filter methylation matrix
colnames_methyl <- colnames(methyl_df)
colnames_methyl_new <- gsub(pattern = "\\-T", replacement = "", x = colnames_methyl)
colnames(methyl_df) <- colnames_methyl_new
methyl_df <- methyl_df[, c("Locus", caseids2test)]
## filter expression matrix
colnames_rna <- colnames(rna_df)
colnames_rna_new <- gsub(pattern = "\\-T", replacement = "", x = colnames_rna)
colnames(rna_df) <- colnames_rna_new
rna_df <- rna_df[rna_df$gene_name %in% genes2test, c("gene_id", "gene_name", caseids2test)]
## see if expression data is duplicated for gene name
rna_df$gene_name[duplicated(rna_df$gene_name)] %>% unique() ## 290 gene names are duplicated
## decide to remove the duplicated gene name because 
## for example those RNA such as Y_RNA, it has 782 gene ids and >2000 probes related, this created a lot more test to do which is is not useful for the correlation with chromatin accessibility
rna_df <- rna_df[!duplicated(rna_df$gene_name),]

# process by gene ---------------------------------------------------------
cor_result_df <- merge(x = probe2gene_df, 
                       y = rna_df %>%
                         select(gene_id, gene_name),
                       by.x = c("gene_HGNC"), by.y = c("gene_name"), all.y = T)
p_vals_vec <- vector(mode = "numeric", length = nrow(cor_result_df))
rhos_vec <- vector(mode = "numeric", length = nrow(cor_result_df))
num_0s.rna_vec <- vector(mode = "numeric", length = nrow(cor_result_df))
num_non_na_vec <- vector(mode = "numeric", length = nrow(cor_result_df))

for (genename_tmp in unique(cor_result_df$gene_HGNC)) {
  idxs_genename <- (cor_result_df$gene_HGNC == genename_tmp)
  ## extract RNA info
  exp_vec <- rna_df[rna_df$gene_name == genename_tmp,caseids2test]; exp_vec <- unlist(exp_vec)
  num_0s <- length(which(exp_vec == 0))
  num_0s.rna_vec[idxs_genename] <- num_0s
  if (num_0s > 83) {
    p_vals_vec[idxs_genename] <- NA
    rhos_vec[idxs_genename] <- NA
    num_non_na_vec[idxs_genename] <- NA
  } else {
    probes2test_tmp <- cor_result_df$probeID[cor_result_df$gene_HGNC == genename_tmp]
    for (probeid_tmp in probes2test_tmp) {
      idx_test <- (idxs_genename & (cor_result_df$probeID == probeid_tmp))
      ## extract methylation info
      methyl_vec <- methyl_df[methyl_df$Locus == probeid_tmp, caseids2test]; methyl_vec <- unlist(methyl_vec)
      ## do spearmen correlation
      test_result <- cor.test(x = methyl_vec, y = exp_vec, method = "spearman")
      p_vals_vec[idx_test] <- test_result$p.value
      rhos_vec[idx_test] <- test_result$estimate
      num_non_na_vec[idx_test] <- length(which(!is.na(methyl_vec)))
    }
  }
  print(paste0("Finished ", which(idxs_genename)))
}
cor_result_df$p_val <- p_vals_vec
cor_result_df$rho <- rhos_vec
cor_result_df$num_0s.rna <- num_0s.rna_vec
cor_result_df$num_non_na <- num_non_na_vec

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Methylation_RNA_Correlation.", run_id, ".tsv")
write.table(x = cor_result_df, file = file2write, quote = F, sep = "\t", row.names = F)



  
  
  