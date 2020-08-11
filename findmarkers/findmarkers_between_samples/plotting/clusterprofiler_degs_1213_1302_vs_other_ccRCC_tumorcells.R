# Yige Wu @WashU Aug 2020

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

# input denpendencies -----------------------------------------------------
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
## input degs table
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_between_samples/findmarkers_1213_1302_vs_others_tumorcells/20200805.v1/findallmarkers_wilcox_1213_1302_vs_others..logfcthreshold0.1.minpct0.1.mindiffpct0.1.tsv")

# convert gene symbol to entrez ids ---------------------------------------
## filter degs
deg_df <- deg_df %>%
  filter(p_val_adj < 0.05)
## get genes to convert
genes2convert <- unique(deg_df$gene)
## retrieve entrezgene_id
genesymbol2entrezid_df <- getBM(attributes=c('entrezgene_id', 'hgnc_symbol'), 
                                filters = 'hgnc_symbol', 
                                values = genes2convert, 
                                mart = ensembl)
## add entrez ids to the deg table
deg_df$entrezgene_id <- mapvalues(x = deg_df$gene, from = genesymbol2entrezid_df$hgnc_symbol, to = as.vector(genesymbol2entrezid_df$entrezgene_id))
## filter markers by entrez id
deg_df <- deg_df %>%
  filter(!is.na(entrezgene_id)) %>%
  filter(entrezgene_id != gene)

# run clusterprofiler -----------------------------------------------------
## initiate result 
ovaresult_df <- NULL
ovaresult_list <- list()
deg_tmp_df <- deg_df
## get gene list
de_genes <- deg_tmp_df$avg_logFC
de_genes
names(de_genes) <- deg_tmp_df$entrezgene_id
## sort
de_genes <- sort(de_genes, decreasing = T)
de_genes

# test over-representation analysis and gene set enrichment using wikipathway ------------------------------
ora_up_tmp <- tryCatch(expr = enricher(gene = names(de_genes), TERM2GENE = wpid2gene, TERM2NAME = wpid2name, pvalueCutoff = 0.1),
                       error = function(e) {warning("ORA failed.");return(NULL)})
if (length(ora_up_tmp) > 0 ) {
  ## convert entrez id to gene symbol
  ora_up_tmp <- setReadable(ora_up_tmp, org.Hs.eg.db, keyType = "ENTREZID")
  ora_up_tmp_df <- as.data.frame(ora_up_tmp)
  
  ## make plots
  p <- cnetplot(x = ora_up_tmp, foldChange=de_genes)
  p <- p + scale_color_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0)
  p <- p + labs(color = "Average Expression\nFold Change")
  file2write <- paste0(dir_out, "cnetplot.png")
  png(filename = file2write, width = 1500, height = 1000, res = 150)
  print(p)
  dev.off()
}

# save outputs ------------------------------------------------------------
## save data frames
file2write <- paste0(dir_out, "ORA.", "Result.", run_id, ".tsv")
write.table(x = ovaresult_df, file = file2write, quote = F, sep = "\t", row.names = F)
## save lists
file2write <- paste0(dir_out, "ORA.", "Result.", run_id, ".RDS")
saveRDS(file = file2write, object = ovaresult_list, compress = T)



