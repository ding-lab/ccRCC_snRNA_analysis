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
## input degs
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/findallmarker_wilcox_tumor_vs_pt_on_katmai/20200929.v1/findallmarkers_wilcox_tumorcells_vs_pt.20200929.v1.tsv")

# convert gene symbol to entrez ids ---------------------------------------
deg_df <- deg_df %>%
  dplyr::rename(genesymbol = row_name) %>%
  filter(genesymbol %in% c("PFKP", "GATM", "PVT1", "ERGIC1", "NDRG1", "LINC00887", "PLIN2", "PAM", 
                                  "TNIP1", "PHYKPL", "LRRC41", "EGLN3", "CDK18",
                                  "SEMA5B", "CIT", "ZNF395", "BIRC3", "C16orf74",
                                  "SPIRE1", "FRMD3", "NNMT", "FAM13A", "GIT2",
                                  "SEMA6A", "FOXP2", "COL23A1", "MXI1", "ST6GAL1",
                                  "EFNA5", "HAVCR1", "LUCAT1", "EVA1C", "KLF7",
                                  "PLOD2", "PLSCR1", "CDH6", "DNAH11", "PDK1",
                                  "PKM", "ZNF608", "GRAMD2B", "TNFAIP8", "LINCO1426",
                                  "SLC38A1", "CP", "SLC6A3", "SGPP2", "ENPP3", 
                                  "MYOCOS", "KCTD3", "TSC22D3", "BARX2", "GAS2L3",
                                  "MSC-AS1", "PTGER3", "SLC16A1-AS1", "HILPDA", "STK39",
                                  "EGFR", "PCSK6", "RETREG1", "LINC02471", "TRIM9",
                                  "TMEM232", "SPAG17", "ABI3BP", "NREP", "RNF145", 
                                  "PRELID2", "ABLIM3", "APBB11P", "MGAM", "TSPAN12", 
                                  "PLCB1"))
  
genes2convert <- unique(deg_df$genesymbol)
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
de_genes <- deg_df$avg_logFC
de_genes
names(de_genes) <- deg_df$entrezgene_id
## sort
de_genes <- sort(de_genes, decreasing = T)
de_genes
## store results
file2write <- paste0(dir_out, "Fold_Changes", ".RDS")
saveRDS(object = de_genes, file = file2write, compress = T)

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
}

# test over-representation analysis and gene set enrichment using wikipathway ------------------------------
enricher_out <- tryCatch(expr = enricher_out <- GSEA(geneList = de_genes, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE, pvalueCutoff = 1),
                         error = function(e) {warning("ORA failed.");return(NULL)})

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

