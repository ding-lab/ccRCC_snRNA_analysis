# Yige Wu @WashU May 2021
## source activate ccrcc_snrna

# set up libraries and output directory -----------------------------------
## getting the path to the current script
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
path_this_script <- thisFile()
## set working directory
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
## library additional libaries
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
library(doParallel)
## set run id
no_cores <- 25
# version_tmp <- "maxCores_minus1"
version_tmp <- paste0(no_cores, "Cores")
# version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input methylation
methyl_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Methylation/CPTAC_ccRCC_discovery_tumor_methylation_betavalue_probe_level_v1.0.tsv")
## input methylation annotation
probe2gene_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/methylation/map_probe_to_gene/20210527.v1/Probe2Gene_HGNC.20210527.v1.tsv")
## input gene expression
rna_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/mRNA/CPTAC_ccRCC_discovery_tumor_mRNA_FPKM_UQ_log2_v1.0.tsv")
## input clear-cell RCC classification
histology_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Case_ID/CPTAC_ccRCC_discovery_caseID_v1.0.tsv")

# filter methylation ------------------------------------------------------
## figure out which cases will be tested
caseids2test <- histology_df$CASE_ID[histology_df$Histologic_Type == "Clear cell renal cell carcinoma"]
length(caseids2test) ## 103
## filter methylation matrix
colnames_methyl <- colnames(methyl_df)
colnames_methyl_new <- gsub(pattern = "\\-T", replacement = "", x = colnames_methyl)
colnames(methyl_df) <- colnames_methyl_new
methyl_df <- methyl_df[, c("Locus", caseids2test)]
## filter
methyl_df <- methyl_df[methyl_df$Locus %in% probe2gene_df$probeID,]
methyl_df <- methyl_df[!duplicated(methyl_df$Locus), ]

# filter RNA --------------------------------------------------------------
## figure out how many genes will be tested
genes2test <- rna_df$gene_name[rna_df$gene_name %in% probe2gene_df$gene_HGNC]
length(genes2test) ## 46206
## filter expression matrix
colnames_rna <- colnames(rna_df)
colnames_rna_new <- gsub(pattern = "\\-T", replacement = "", x = colnames_rna)
colnames(rna_df) <- colnames_rna_new
rna_df <- rna_df[rna_df$gene_name %in% genes2test, c("gene_id", "gene_name", caseids2test)]
## see if expression data is duplicated for gene name
# rna_df$gene_name[duplicated(rna_df$gene_name)] %>% unique() ## 290 gene names are duplicated
## decide to remove the duplicated gene name because 
## for example those RNA such as Y_RNA, it has 782 gene ids and >2000 probes related, this created a lot more test to do which is is not useful for the correlation with chromatin accessibility
rna_df <- rna_df[!duplicated(rna_df$gene_name),]
genes2test <- rna_df$gene_nam
length(genes2test) ## 44391

# process by gene ---------------------------------------------------------
registerDoParallel(cores = no_cores)
cor_result_list<-foreach(i=1:length(genes2test)) %dopar% {
  ## extract RNA info
  g <- genes2test[i]
  exp_vec <- rna_df[rna_df$gene_name == g,3:ncol(rna_df)]; exp_vec <- unlist(exp_vec)
  num_0s_tmp <- length(which(exp_vec == 0))
  probes2test_tmp <- probe2gene_df$probeID[probe2gene_df$gene_HGNC == g]
  test_result_df <- NULL
  for (probeid_tmp in probes2test_tmp) {
    ## extract methylation info
    methyl_vec <- methyl_df[methyl_df$Locus == probeid_tmp,2:ncol(methyl_df)]; methyl_vec <- unlist(methyl_vec)
    ## do spearmen correlation
    test_result <- cor.test(x = methyl_vec, y = exp_vec, method = "spearman")
    test_result_tmp_df <- data.frame(probeID = probeid_tmp, p_val = test_result$p.value, rho = test_result$estimate, num_non_na = length(which(!is.na(methyl_vec))), num_0s = num_0s_tmp)
    test_result_df <- rbind(test_result_df, test_result_tmp_df)
  }
  test_result_df$gene_symbol = g
  print(paste0(i, ":", signif(x = 100*(i/length(genes2test)), digits = 4), "% done"))
  return(test_result_df)
}
file2write <- paste0(dir_out, "Methylation_RNA_Correlation.", run_id, ".RDS")
saveRDS(object = cor_result_list, file = file2write, compress = T)
cor_result_df <- do.call(rbind.data.frame, cor_result_list)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Methylation_RNA_Correlation.", run_id, ".tsv")
write.table(x = cor_result_df, file = file2write, quote = F, sep = "\t", row.names = F)


