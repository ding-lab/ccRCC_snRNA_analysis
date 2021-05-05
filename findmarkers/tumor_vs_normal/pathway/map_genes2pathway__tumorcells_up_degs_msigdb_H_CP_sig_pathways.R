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
## input the genes to plot
pathway2genes_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/pathway/ora_msigdb_H_CP_tumorcells_up_degs/20210429.v2/ORA_Results.tsv")

# get gene-pathway map-------------------------------------------------
## add name for the marker groups
pathway2genes_filtered_df <- pathway2genes_df %>%
  filter(p.adjust < 0.05)
pathway2genes_list <- sapply(pathway2genes_filtered_df$geneID, function(x) {
  genes_vec <- str_split(string = x, pattern = "\\/")[[1]]
  return(genes_vec)
})
gene2pathway_df <- data.frame(GeneSet_Name = rep(x = pathway2genes_filtered_df$Description, sapply(X = pathway2genes_list, FUN = function(x) {
  length_vec <- length(x)
  return(length_vec)
})), GeneSymbol = unlist(pathway2genes_list))
gene2pathway_df <- unique(gene2pathway_df)
gene2pathways_df <- gene2pathway_df %>%
  group_by(GeneSymbol) %>%
  summarise(GeneSet_Names = paste0(sort(GeneSet_Name), collapse = "|")) %>%
  arrange(GeneSet_Names)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "DEG2Pathway.tsv")
write.table(x = gene2pathway_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "DEG2Pathway.Collapsed.tsv")
write.table(x = gene2pathways_df, file = file2write, quote = F, sep = "\t", row.names = F)

