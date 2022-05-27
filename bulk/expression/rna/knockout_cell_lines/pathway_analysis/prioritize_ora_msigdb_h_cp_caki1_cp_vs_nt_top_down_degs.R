# Yige Wu @ WashU 2022 May

# set up libraries and output directory -----------------------------------
## set working directory
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_Drug/"
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    install.packages(pkg_name_tmp, dependencies = T)
  }
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    if (!requireNamespace("BiocManager", quietly=TRUE))
      install.packages("BiocManager")
    BiocManager::install(pkg_name_tmp)
  }
  library(package = pkg_name_tmp, character.only = T)
}

# input -------------------------------------------------------------------
enricher_all_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/rna/knockout_cell_lines/pathway_analysis/ora_msigdb_h_cp_caki1_cp_vs_nt_top_down_degs/20220517.v1/ORA_Results.tsv")
enricher_all_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/rna/knockout_cell_lines/pathway_analysis/ora_msigdb_h_cp_caki1_cp_vs_nt_top_down_degs/20220517.v2/ORA_Results.tsv")

# preprocess -----------------------------------------------------------------
enricher_all_df <- enricher_all_df %>%
  mutate(number_genes_tested = str_split_fixed(string = GeneRatio, pattern = "\\/", n = 2)[,2]) %>%
  mutate(number_genes_tested = as.numeric(number_genes_tested)) %>%
  mutate(generatio_num = Count/number_genes_tested)%>%
  mutate(row_id = ID) %>%
  mutate(test = "Caki1_cp_vs_nt") ## need to add this: test group
enricher_filtered_df <- enricher_all_df %>%
  filter(p.adjust < 0.05) %>%
  arrange(desc(generatio_num))
enricher_filtered_df$Keep <- F

# process by manual repitation ---------------------------------------------------------
## step 1: add one pathway at a time based on the max overlap
### for each group, add the top pathway first and add one for each iteration
#### top 500 genes
# enricher_filtered_df$Keep[enricher_filtered_df$ID %in% c("NABA_MATRISOME_ASSOCIATED", "REACTOME_SIGNALING_BY_INTERLEUKINS", 
#                                                          "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_KRAS_SIGNALING_UP",
#                                                          "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "REACTOME_BIOLOGICAL_OXIDATIONS",
#                                                          "HALLMARK_COMPLEMENT", "REACTOME_FORMATION_OF_THE_CORNIFIED_ENVELOPE",
#                                                          "WP_SELENIUM_MICRONUTRIENT_NETWORK", "REACTOME_INTERFERON_GAMMA_SIGNALING")] <- T
#### top 1000 genes
# enricher_filtered_df$Keep[enricher_filtered_df$ID %in% c("NABA_MATRISOME_ASSOCIATED", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
#                                                          "HALLMARK_INTERFERON_GAMMA_RESPONSE", "REACTOME_SIGNALING_BY_INTERLEUKINS",
#                                                          "REACTOME_HEMOSTASIS", "HALLMARK_INFLAMMATORY_RESPONSE",
#                                                          "HALLMARK_KRAS_SIGNALING_UP", "REACTOME_SIGNALING_BY_GPCR",
#                                                          "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "NABA_CORE_MATRISOME")] <- T
enricher_filtered_df$Keep[enricher_filtered_df$ID %in% c("NABA_MATRISOME_ASSOCIATED", "HALLMARK_INTERFERON_GAMMA_RESPONSE",
                                                         "REACTOME_HEMOSTASIS", "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
                                                         "REACTOME_SIGNALING_BY_INTERLEUKINS", "HALLMARK_INFLAMMATORY_RESPONSE",
                                                         "REACTOME_SIGNALING_BY_GPCR", "HALLMARK_KRAS_SIGNALING_UP",
                                                         "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "REACTOME_NEUTROPHIL_DEGRANULATION")] <- T
## step 2: run this
enricher_filtered_df$max_overlap_ratio <- sapply(enricher_filtered_df$row_id, function(row_id_tmp, test_df) {
  # sapply(head(enricher_filtered_df$row_id), function(row_id_tmp, test_df) {
  
  test_tmp <- test_df$test[test_df$row_id == row_id_tmp]
  generatio_tmp <- test_df$generatio_num[test_df$row_id == row_id_tmp]
  id_tmp <- test_df$ID[test_df$row_id == row_id_tmp]
  test_keep_df <- test_df[test_df$test == test_tmp & test_df$Keep & test_df$generatio_num >= generatio_tmp & test_df$ID != id_tmp,]
  
  genes_tmp <- str_split(string = test_df$geneID[test_df$row_id == row_id_tmp], pattern = "\\/")[[1]]
  max_overlap_ratio <- 0
  for (genestring_kept_tmp in test_keep_df$geneID) {
    genes_kept_tmp <- str_split(string = genestring_kept_tmp, pattern = "\\/")[[1]]
    genes_common_tmp <- intersect(genes_tmp, genes_kept_tmp)
    overlap_ratio_tmp <- length(genes_common_tmp)/length(genes_tmp)
    max_overlap_ratio <- max(c(max_overlap_ratio, overlap_ratio_tmp))
  }
  return(max_overlap_ratio)
}, test_df = enricher_filtered_df)
## scroll down the table and select the ones with max overlap % < 50% and meaningful
View(enricher_filtered_df)

# output ------------------------------------------------------------------
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis//functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

file2write <- paste0(dir_out, "ora_msigdb_h_cp_caki1_cp_vs_nt_top_down_degs.processed.", run_id, ".tsv")
write.table(x = enricher_all_df, file = file2write, quote = F, sep = "\t", row.names = F)

file2write <- paste0(dir_out, "ora_msigdb_h_cp_caki1_cp_vs_nt_top_down_degs.processed.sig.", run_id, ".tsv")
write.table(x = enricher_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)

