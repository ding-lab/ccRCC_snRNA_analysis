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
pathway2genes_df1 <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/bap1/pathway/ora_msigdb_H_CP_BAP1_vs_NonMutant_down_daps_promoter_enhancer_28samples/20211011.v1/ORA_Results.tsv")
pathway2genes_df2 <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/bap1/pathway/ora_msigdb_H_CP_BAP1_vs_NonMutant_up_daps_promoter_enhancer_28samples/20211011.v1/ORA_Results.tsv")
pathway2genes_df3 <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/pbrm1/pathway/ora_msigdb_H_CP_PBRM1_vs_NonMutant_down_daps_promoter_enhancer_28samples/20211011.v1/ORA_Results.tsv")
pathway2genes_df4 <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/pbrm1/pathway/ora_msigdb_H_CP_PBRM1_vs_NonMutant_up_daps_promoter_enhancer_28samples/20211011.v1/ORA_Results.tsv")
## input pathway pairwise comparison
# ora_obj1 <- readRDS(file = "./Resources/Analysis_Results/snatac/da_peaks/bap1/pathway/ora_msigdb_H_CP_BAP1_vs_NonMutant_down_daps_promoter_enhancer/20210625.v1/ORA_Results.RDS")
# ora_obj4 <- readRDS(file = "./Resources/Analysis_Results/snatac/da_peaks/bap1/pathway/ora_msigdb_H_CP_PBRM1_vs_NonMutant_up_daps_promoter_enhancer/20210625.v1/ORA_Results.RDS")

# preprocess pairwise comparison ------------------------------------------
similarity_list <- list()
ora_obj4 <- enrichplot::pairwise_termsim(ora_obj4)
similarity_mat4 <- ora_obj4@termsim
rm(ora_obj4)
ora_obj1 <- enrichplot::pairwise_termsim(ora_obj1)
similarity_mat1 <- ora_obj1@termsim
rm(ora_obj1)
similarity_list[["PBRM1_Up"]] <- similarity_mat4
similarity_list[["BAP1_Down"]] <- similarity_mat1

# unite and select pathways -----------------------------------------------------------------
pathway2genes_df <- rbind(pathway2genes_df1 %>%
                            mutate(Comparison = "BAP1_Down"),
                          pathway2genes_df2 %>%
                            mutate(Comparison = "BAP1_Up"),
                          pathway2genes_df3 %>%
                            mutate(Comparison = "PBRM1_Down"),
                          pathway2genes_df4 %>%
                            mutate(Comparison = "PBRM1_Up"))
pathway2genes_sig_df <- pathway2genes_df[pathway2genes_df$p.adjust < 0.05,]
## select pathways
### select top 3 most significantly enriched pathways (jaccard similarity coefficient < 0.1) for each comparisons
pathway_selected <- c("HALLMARK_MTORC1_SIGNALING", "WP_EGFEGFR_SIGNALING_PATHWAY","REACTOME_RHO_GTPASE_CYCLE", "WP_DNA_REPAIR_PATHWAYS_FULL_NETWORK",
                      "WP_FOCAL_ADHESION", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "PID_EPHA_FWDPATHWAY", "WP_TGFB_SIGNALING_IN_THYROID_CELLS_FOR_EPITHELIALMESENCHYMAL_TRANSITION",
                      "WP_FOCAL_ADHESION", "REACTOME_RAC1_GTPASE_CYCLE", "HALLMARK_HYPOXIA")

pathway2genes_filtered_df <- pathway2genes_df %>%
  filter(Description %in% pathway_selected) %>%
  group_by(Comparison) %>%
  mutate(rank_pvalue = rank(pvalue, ties.method = "first"))

# map each gene to pathway ------------------------------------------------
pathway2genes_list <- sapply(pathway2genes_filtered_df$geneID, function(x) {
  genes_vec <- str_split(string = x, pattern = "\\/")[[1]]
  return(genes_vec)
})
genes2filter <- unique(unlist(pathway2genes_list))
gene2pathway_df <- data.frame(GeneSet_Name = rep(x = pathway2genes_filtered_df$Description, sapply(X = pathway2genes_list, FUN = function(x) {
  length_vec <- length(x)
  return(length_vec)
})), GeneSymbol = unlist(pathway2genes_list))
gene2pathway_df <- unique(gene2pathway_df)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "PBRM1_BAP1_vs_NonMutants.DAPGene2TopPathway.", run_id, ".tsv")
write.table(x = gene2pathway_df, file = file2write, sep = "\t", row.names = F, quote = F)
file2write <- paste0(dir_out, "PBRM1_BAP1_vs_NonMutants.DAPPathways2Genes.", run_id, ".tsv")
write.table(x = pathway2genes_df, file = file2write, sep = "\t", row.names = F, quote = F)
file2write <- paste0(dir_out, "PBRM1_BAP1_vs_NonMutants.DAPTopPathways2Genes.", run_id, ".tsv")
write.table(x = pathway2genes_filtered_df, file = file2write, sep = "\t", row.names = F, quote = F)
# file2write <- paste0(dir_out, "PBRM1_BAP1_vs_NonMutants.DAPPathwaysSimilarity.", run_id, ".RDS")
# saveRDS(object = similarity_list, file = file2write, compress = T)
