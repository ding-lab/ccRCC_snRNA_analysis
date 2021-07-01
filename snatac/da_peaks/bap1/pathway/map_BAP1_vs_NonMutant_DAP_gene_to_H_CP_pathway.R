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
pathway2genes_df1 <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/bap1/pathway/ora_msigdb_H_CP_BAP1_vs_NonMutant_down_daps_promoter_enhancer/20210628.v1/ORA_Results.tsv")
pathway2genes_df2 <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/bap1/pathway/ora_msigdb_H_CP_BAP1_vs_NonMutant_up_daps_promoter_enhancer/20210628.v1/ORA_Results.tsv")
## input meta data
id_metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210322.v1/meta_data.20210322.v1.tsv")
## input cell number per cluster
genesetnames_plot <- c("HALLMARK_MTORC1_SIGNALING", "WP_EGFEGFR_SIGNALING_PATHWAY","REACTOME_RHO_GTPASE_CYCLE", "WP_DNA_REPAIR_PATHWAYS_FULL_NETWORK",
                      "WP_FOCAL_ADHESION", "HALLMARK_TNFA_SIGNALING_VIA_NFKB", "PID_EPHA_FWDPATHWAY", "WP_TGFB_SIGNALING_IN_THYROID_CELLS_FOR_EPITHELIALMESENCHYMAL_TRANSITION",
                      "WP_FOCAL_ADHESION", "REACTOME_RAC1_GTPASE_CYCLE", "HALLMARK_HYPOXIA")
genesetnames_plot <- unique(genesetnames_plot)


# preprocess --------------------------------------------------------------
pathway2genes_df <- rbind(pathway2genes_df1, pathway2genes_df2)

# get gene-pathway map-------------------------------------------------
## add name for the marker groups
pathway2genes_filtered_df <- pathway2genes_df %>%
  filter(Description %in% genesetnames_plot)
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
file2write <- paste0(dir_out, "BAP1_vs_NonMutant.DAP.Gene2Pathway.tsv")
write.table(x = gene2pathway_df, file = file2write, quote = F, sep = "\t", row.names = F)
