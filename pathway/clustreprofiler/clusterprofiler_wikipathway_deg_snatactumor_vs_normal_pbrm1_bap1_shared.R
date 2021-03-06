# Yige Wu @WashU Oct 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
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
## input the deg annotation
deg_anno_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/annotate_deg/annotate_deg_by_snatactumorgroup_shared/20201130.v1/DEG_each_tumor_vs_pt.annotated.tsv")
## get biomart
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
## input wikipathway 
wp2gene <- read.gmt("./Resources/Knowledge/Databases/Wikipathways/wikipathways-20200510-gmt-Homo_sapiens.gmt")
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

# convert gene symbol to entrez ids ---------------------------------------
deg_df <- deg_anno_df %>%
  filter(category_byshared %in% c("FALSE_TRUE_TRUE"))
genes2convert <- unique(deg_df$genesymbol_deg)
## retrieve entrezgene_id
genesymbol2entrezid_df <- getBM(attributes=c('entrezgene_id', 'hgnc_symbol'), 
                                filters = 'hgnc_symbol', 
                                values = genes2convert, 
                                mart = ensembl)
## add entrez ids to the deg table
deg_df$entrezgene_id <- mapvalues(x = deg_df$genesymbol, from = genesymbol2entrezid_df$hgnc_symbol, to = as.vector(genesymbol2entrezid_df$entrezgene_id))
## filter markers by entrez id
deg_df <- deg_df %>%
  filter(!is.na(entrezgene_id)) %>%
  filter(entrezgene_id != genesymbol)

# prepare inputs for clusterprofiler -----------------------------------------------------
## get gene list
de_genes <- deg_df$mean_avg_logFC.bap1mutant
de_genes
names(de_genes) <- deg_df$entrezgene_id
## sort
de_genes <- sort(de_genes, decreasing = T)
de_genes
## store results
file2write <- paste0(dir_out, "Fold_Changes", ".RDS")
saveRDS(object = de_genes, file = file2write, compress = T)

# test over-representation analysis and gene set enrichment using wikipathway ------------------------------
enricher_out <- tryCatch(expr = enricher(gene = names(de_genes), TERM2GENE = wpid2gene, TERM2NAME = wpid2name, qvalueCutoff = 1, pvalueCutoff = 1),
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
}

# test over-representation analysis and gene set enrichment using wikipathway ------------------------------
enricher_out <- tryCatch(expr = enricher_out <- GSEA(geneList = de_genes, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, qvalueCutoff = 1, pvalueCutoff = 1),
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
}

