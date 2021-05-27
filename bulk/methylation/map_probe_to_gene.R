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

# filter methylation ------------------------------------------------------
nrow(methyl_df) ## 917437 rows
## some of probes are duplicated, and the values are the same across samples
methyl_df$Locus[!(methyl_df$Locus %in% methyl_anno_df$probeID)]
methyl_df$Locus[duplicated(methyl_df$Locus)]
# methyl_df[methyl_df$Locus == "cg00000714",]
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
