# Yige Wu @WashU May 2020

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
markers_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_with_cna/findallmarkers_genelevel_expected_cnv_vs_neutral_in_tumorcells/20200518.v1/ExpectedCNV_vs_Neutral..FindAllMarkers.Wilcox..20200518.v1.tsv")

# convert gene symbol to entrez ids ---------------------------------------
## get biomart
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
## filter markers
markers_df <- markers_df %>%
  filter(p_val_adj < 0.05)
## retrieve UCSC stable ids
genesymbol2ucsc_df <- getBM(attributes=c('entrezgene_id', 'hgnc_symbol'), 
                            filters = 'hgnc_symbol', 
                            values = unique(markers_df$de_gene_symbol), 
                            mart = ensembl)
## add entrez ids to the deg table
markers_df$entrezgene_id <- mapvalues(x = markers_df$de_gene_symbol, from = genesymbol2ucsc_df$hgnc_symbol, to = as.vector(genesymbol2ucsc_df$entrezgene_id))
## filter markers by entrez id
markers_df <- markers_df %>%
  filter(!is.na(entrezgene_id)) %>%
  filter(entrezgene_id != de_gene_symbol)

# process gmt -------------------------------------------------------------
wp2gene <- read.gmt("./Resources/Databases/Wikipathways/wikipathways-20200510-gmt-Homo_sapiens.gmt")
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

# make ranked gene list ---------------------------------------------------
cna_gene_plot <- "VHL"
## get genes
markers_for_cnagene_df <- markers_df %>%
  filter(cna_gene_symbol == cna_gene_plot) %>%
  group_by(de_gene_symbol, entrezgene_id) %>%
  summarise(num_aliquot_de = n(), mean_avg_logFC = mean(avg_logFC))
markers_for_cnagene_df <- as.data.frame(markers_for_cnagene_df)
## get gene list
de_genes <- markers_for_cnagene_df$mean_avg_logFC
de_genes
names(de_genes) <- markers_for_cnagene_df$entrezgene_id
## sort
de_genes <- sort(de_genes, decreasing = T)
de_genes

# test over-representation analysis and gene set enrichment using wikipathway ------------------------------
ewp <- enricher(gene = names(de_genes), TERM2GENE = wpid2gene, TERM2NAME = wpid2name)
if (nrow(ewp) > 0) {
  ewp <- setReadable(ewp, org.Hs.eg.db, keyType = "ENTREZID")
  ewp_df <- as.data.frame(ewp)
}
ewp2 <- GSEA(geneList = de_genes, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE, pvalueCutoff = 0.1)
if (nrow(ewp2) > 0) {
  ewp2 <- setReadable(ewp2, org.Hs.eg.db, keyType = "ENTREZID")
  ewp2_df <- as.data.frame(ewp2)
}
