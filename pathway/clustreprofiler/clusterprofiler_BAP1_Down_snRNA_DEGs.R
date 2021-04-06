# Yige Wu @WashU Sep 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(clusterProfiler)
library(biomaRt)
library(org.Hs.eg.db)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## get biomart
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
## input wikipathway 
wp2gene <- read.gmt("./Resources/Knowledge/Databases/Wikipathways/wikipathways-20200510-gmt-Homo_sapiens.gmt")
wp2gene <- wp2gene %>% tidyr::separate(term, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME
## input degs
### 2021-04-05
deg_all_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/pbrm1_bap1_vs_non_mutants/summarize_degs/summarize_PBRM1_BAP1_DEGs/20210405.v1/BAP1_DEGs.20210405.v1.tsv")
avg_logFC_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/pbrm1_bap1_vs_non_mutants/summarize_degs/summarize_PBRM1_BAP1_DEGs/20210405.v1/BAP1_PBRM1_DEGs.Mean_avg_logFC.20210405.v1.tsv")

# convert gene symbol to entrez ids ---------------------------------------
deg_df <- deg_all_df %>%
  filter(BAP1_vs_NonMutants_snRNA == "Down")
genes2convert <- unique(deg_df$genesymbol_deg)
## retrieve entrezgene_id
genesymbol2entrezid_df <- getBM(attributes=c('entrezgene_id', 'hgnc_symbol'), 
                                filters = 'hgnc_symbol', 
                                values = genes2convert, 
                                mart = ensembl)
## add entrez ids to the deg table
deg_df$entrezgene_id <- mapvalues(x = deg_df$genesymbol_deg, from = genesymbol2entrezid_df$hgnc_symbol, to = as.vector(genesymbol2entrezid_df$entrezgene_id))
## filter markers by entrez id
deg_df <- deg_df %>%
  filter(!is.na(entrezgene_id)) %>%
  filter(entrezgene_id != genesymbol_deg)

# prepare inputs for clusterprofiler -----------------------------------------------------
## get gene list
avg_logFC_df$avg_logFC_BAP1_vs_others <- rowMeans(avg_logFC_df[, c("BAP1_Tumorcells_vs_PTcells_Down", "BAP1_vs_NonMutants_Tumorcells_Down", "BAP1_vs_PBRM1_Mutants_Tumorcells_Down")])
deg_df$avg_logFC <- mapvalues(x = deg_df$genesymbol_deg, from = avg_logFC_df$genesymbol_deg, to = as.vector(avg_logFC_df$avg_logFC_BAP1_vs_others))
deg_df$avg_logFC <- as.numeric(deg_df$avg_logFC)
de_genes <- -(deg_df$avg_logFC)
de_genes
names(de_genes) <- deg_df$entrezgene_id
## sort
de_genes <- sort(de_genes, decreasing = T)
de_genes

# test over-representation analysis and gene set enrichment using wikipathway ------------------------------
enricher_out <- tryCatch(expr = enricher(gene = names(de_genes), TERM2GENE = wpid2gene, TERM2NAME = wpid2name, pvalueCutoff = 1),
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
  png(filename = file2write, width = 2000, height = 1500, res = 150)
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


# save output -------------------------------------------------------------
## store results
file2write <- paste0(dir_out, "Fold_Changes", ".RDS")
saveRDS(object = de_genes, file = file2write, compress = T)
