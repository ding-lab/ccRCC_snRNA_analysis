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
## input probe annotation file
probe_anno_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Methylation/EPIC.hg38.manifest.tsv")
## input the original methylation data
methyl_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Methylation/CCRCC.TumorBeta.txt")
## input BAP1 DACRs
bap1_dacr_df <- fread(data.table = F, input = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/Differential_Peaks/BAP1_specific_DACRs_Promoter.20210526.tsv")

# process --------------------------------------------------------------
## filter
probe_anno_df <- probe_anno_df %>%
  filter(!is.na(gene))

## split by gene
idx_rep <- sapply(probe_anno_filtered_df$gene, function(x) {
  genes_vec <- str_split(string = x, pattern = ";")[[1]]
  genes_len <- length(genes_vec)
  return(genes_len)
})
probe2gene_df <- probe_anno_filtered_df[rep(1:nrow(probe_anno_filtered_df), idx_rep),]
probe2gene_df$gene <- unlist(sapply(probe_anno_filtered_df$gene, function(x) {
  genes_vec <- str_split(string = x, pattern = ";")[[1]]
  return(genes_vec)
}))
probe2gene_df$gene_HGNC <- unlist(sapply(probe_anno_filtered_df$gene_HGNC, function(x) {
  genes_vec <- str_split(string = x, pattern = ";")[[1]]
  return(genes_vec)
}))
  
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "methyl_subtype_specific_1000_probe2allgenes.", run_id, ".tsv")
write.table(x = probe_anno_filtered_df, file = file2write, sep = "\t", row.names = F, quote = F)
file2write <- paste0(dir_out, "methyl_subtype_specific_1000_probe2gene.", run_id, ".tsv")
write.table(x = probe2gene_df, file = file2write, sep = "\t", row.names = F, quote = F)
