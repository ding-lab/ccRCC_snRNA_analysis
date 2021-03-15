# Yige Wu @WashU Sep 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
if (!requireNamespace("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
  BiocManager::install("org.Hs.eg.db")
if (!requireNamespace("magrittr", quietly = TRUE))
  install.packages("magrittr")
library(magrittr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(ggnewscale)
library(ggplot2)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input wikipathway 
wp2gene <- clusterProfiler::read.gmt("./Resources/Knowledge/Databases/Wikipathways/wikipathways-20200510-gmt-Homo_sapiens.gmt")
wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
## input degs
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/edgeR/run_deg_analysis/run_PBRM1_BAP1_deg_on_cptac_ccRCC_discovery_cases/20210312.v2/PBRM1_BAP1_DEGs.glmQLFTest.outputtables.tsv")
## input the gene id mapping
gene_mapping_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/edgeR/map_to_genesymbols/map_CPTAC_Discovery_ensemblgeneids_to_genesymbols/20210312.v1/ensembl_gene_id.mapping.nodup.20210312.v1.tsv")

# preprocess ----------------------------------------------------------
## filter
deg_filtered_df <- deg_df %>%
  filter(logFC > 0) %>%
  filter(FDR < 0.05) %>%
  filter(grepl(pattern = "vs Non-mutants", x = comparison)) %>%
  filter(!(gene_ensembl_id %in% c("__alignment_not_unique", "__ambiguous", "__no_feature")))
logFC_wide_df <- dcast(data = deg_filtered_df, formula = gene_ensembl_id ~ comparison, value.var = "logFC")
genes2convert_df <- logFC_wide_df %>%
  filter(!is.na(`BAP1 mutated vs Non-mutants`) & is.na(`PBRM1 mutated vs Non-mutants`))
  

# convert gene symbol to entrez ids ---------------------------------------
genes2convert_df <- genes2convert_df %>%
  dplyr::rename(ensembl_gene_id_version.deg = gene_ensembl_id) %>%
  mutate(ensembl_gene_id = str_split_fixed(string = ensembl_gene_id_version.deg, pattern = "\\.", n = 2)[,1])
genes2convert_df <- merge(x = genes2convert_df, y = gene_mapping_df, by = c("ensembl_gene_id"), all.x = T)

# prepare inputs for clusterprofiler -----------------------------------------------------
genes2convert_filtered_df <- genes2convert_df %>%
  filter(!is.na(entrezgene_id))
## get gene list
de_genes <- genes2convert_filtered_df$`BAP1 mutated vs Non-mutants`
de_genes
names(de_genes) <- genes2convert_filtered_df$entrezgene_id
## sort
de_genes <- sort(de_genes, decreasing = T)
de_genes
## store results
file2write <- paste0(dir_out, "Fold_Changes", ".RDS")
saveRDS(object = de_genes, file = file2write, compress = T)
## make universe
entrezids_universe <- gene_mapping_df$entrezgene_id[!is.na(gene_mapping_df$entrezgene_id)]
entrezids_universe <- as.character(entrezids_universe)

# test over-representation analysis and gene set enrichment using wikipathway ------------------------------
enricher_out <- tryCatch(expr = enricher(gene = names(de_genes), TERM2GENE = wpid2gene, TERM2NAME = wpid2name, pvalueCutoff = 1, universe = entrezids_universe),
                         error = function(e) {warning("ORA failed.");return(NULL)})

if (length(enricher_out) > 0 ) {
  ## convert entrez id to gene symbol
  enricher_out <- setReadable(enricher_out, org.Hs.eg.db, keyType = "ENTREZID")
  enricher_out_df <- as.data.frame(enricher_out)
  ## store results
  file2write <- paste0(dir_out, "ORA_Results", ".RDS")
  saveRDS(object = enricher_out, file = file2write, compress = T)
  ## store results
  file2write <- paste0(dir_out, "ORA_Results", ".tsv")
  write.table(x = enricher_out_df, file = file2write, quote = F, row.names = F, sep = "\t")
  ## make plots
  p <- cnetplot(x = enricher_out, foldChange=de_genes)
  p <- p + scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0)
  p <- p + labs(color = "Average Expression\nFold Change")
  file2write <- paste0(dir_out, "ORA", ".png")
  png(filename = file2write, width = 1200, height = 1000, res = 150)
  print(p)
  dev.off()
}

# save output -------------------------------------------------------------
# test over-representation analysis and gene set enrichment using wikipathway ------------------------------
enricher_out <- tryCatch(expr = enricher_out <- GSEA(geneList = de_genes, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE, pvalueCutoff = 1),
                         error = function(e) {warning("GSEA failed.");return(NULL)})

if (length(enricher_out) > 0 ) {
  ## convert entrez id to gene symbol
  enricher_out <- setReadable(enricher_out, org.Hs.eg.db, keyType = "ENTREZID")
  enricher_out_df <- as.data.frame(enricher_out)
  
  ## store results
  file2write <- paste0(dir_out, "GSEA_Results", ".RDS")
  saveRDS(object = enricher_out, file = file2write, compress = T)
  ## store results
  file2write <- paste0(dir_out, "GSEA_Results", ".tsv")
  write.table(x = enricher_out_df, file = file2write, quote = F, row.names = F, sep = "\t")
  ## make plots
  p <- cnetplot(x = enricher_out, foldChange=de_genes)
  p <- p + scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0)
  p <- p + labs(color = "Average Expression\nFold Change")
  file2write <- paste0(dir_out, "GSEA", ".png")
  png(filename = file2write, width = 1200, height = 1000, res = 150)
  print(p)
  dev.off()
}

