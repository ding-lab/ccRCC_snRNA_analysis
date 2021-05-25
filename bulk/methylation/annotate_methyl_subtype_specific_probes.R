# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input methylylation probe group
probe_group_df1 <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Methylation/Methylation_Subtype/methylation.sig.50.probes.by.metSubtype_3.1.txt")
probe_group_df2 <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Methylation/Methylation_Subtype/methylation.sig.50.probes.by.metSubtype_3.2.txt")
probe_group_df3 <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Methylation/Methylation_Subtype/methylation.sig.50.probes.by.metSubtype_3.3.txt")
## input probe annotation file
probe_anno_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Methylation/EPIC.hg38.manifest.tsv")

# process --------------------------------------------------------------
probe_group_df <- rbind(probe_group_df1 %>%
                          mutate(Probe_subtype = "m1"),
                        probe_group_df2 %>%
                          mutate(Probe_subtype = "m2"))
probe_group_df <- rbind(probe_group_df,
                        probe_group_df3 %>%
                          mutate(Probe_subtype = "m3"))
colnames(probe_group_df)[1] <- "probeID"
## filter
probe_anno_df <- probe_anno_df %>%
  filter(probeID %in% probe_group_df$probeID)
## merge
probe_group_df$gene <- mapvalues(x = probe_group_df$probeID, from = probe_anno_df$probeID, to = as.vector(probe_anno_df$gene))
probe_group_df$gene_HGNC <- mapvalues(x = probe_group_df$probeID, from = probe_anno_df$probeID, to = as.vector(probe_anno_df$gene_HGNC))
probe_group_df$CpG_chrm <- mapvalues(x = probe_group_df$probeID, from = probe_anno_df$probeID, to = as.vector(probe_anno_df$CpG_chrm))
probe_group_df$probeBeg <- mapvalues(x = probe_group_df$probeID, from = probe_anno_df$probeID, to = as.vector(probe_anno_df$probeBeg))
probe_group_df$probeEnd <- mapvalues(x = probe_group_df$probeID, from = probe_anno_df$probeID, to = as.vector(probe_anno_df$probeEnd))
probe_group_df$CpG_beg <- mapvalues(x = probe_group_df$probeID, from = probe_anno_df$probeID, to = as.vector(probe_anno_df$CpG_beg))
probe_group_df$CpG_end <- mapvalues(x = probe_group_df$probeID, from = probe_anno_df$probeID, to = as.vector(probe_anno_df$CpG_end))
## split by gene
idx_rep <- sapply(probe_group_df$gene, function(x) {
  genes_vec <- str_split(string = x, pattern = ";")[[1]]
  genes_len <- length(genes_vec)
  return(genes_len)
})
probe2gene_df <- probe_group_df[rep(1:nrow(probe_group_df), idx_rep),]
probe2gene_df$gene <- unlist(sapply(probe_group_df$gene, function(x) {
  genes_vec <- str_split(string = x, pattern = ";")[[1]]
  return(genes_vec)
}))
probe2gene_df$gene_HGNC <- unlist(sapply(probe_group_df$gene_HGNC, function(x) {
  genes_vec <- str_split(string = x, pattern = ";")[[1]]
  return(genes_vec)
}))
  
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "methyl_subtype_specific_probe2allgenes.", run_id, ".tsv")
write.table(x = probe_group_df, file = file2write, sep = "\t", row.names = F, quote = F)
file2write <- paste0(dir_out, "methyl_subtype_specific_probe2gene.", run_id, ".tsv")
write.table(x = probe2gene_df, file = file2write, sep = "\t", row.names = F, quote = F)
