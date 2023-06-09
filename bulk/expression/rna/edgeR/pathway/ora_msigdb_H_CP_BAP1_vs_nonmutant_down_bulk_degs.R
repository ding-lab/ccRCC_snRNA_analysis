# Yige Wu @WashU Sep 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(clusterProfiler)
library(org.Hs.eg.db)
library(biomaRt)

## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input wikipathway 
wp2gene1 <- read.gmt("./Resources/Knowledge/Databases/MSigDB/h.all.v7.4.entrez.gmt")
wp2gene2 <- read.gmt("./Resources/Knowledge/Databases/MSigDB/c2.cp.v7.4.entrez.gmt")
wp2gene <- rbind(wp2gene1, wp2gene2)
## input degs
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/edgeR/run_deg_analysis/run_PBRM1_BAP1_deg_on_cptac_ccRCC_discovery_cases/20210329.v1/PBRM1_BAP1_DEGs.glmQLFTest.OutputTables.tsv")

# convert gene symbol to entrez ids ---------------------------------------
genes_background_df <- deg_df
genes2convert <- unique(genes_background_df$ensembl_gene_id)
## retrieve entrezgene_id
ensembl <- useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)

genesymbol2entrezid_df <- getBM(attributes=c('entrezgene_id', 'hgnc_symbol', "ensembl_gene_id"), 
                                filters = 'ensembl_gene_id', 
                                values = genes2convert, 
                                mart = ensembl)
## add entrez ids to the deg table
genes_background_df$entrezgene_id <- mapvalues(x = genes_background_df$ensembl_gene_id, from = genesymbol2entrezid_df$ensembl_gene_id, to = as.vector(genesymbol2entrezid_df$entrezgene_id))
##
genes_test_df <- genes_background_df %>%
  filter(logFC < -1 & FDR < 0.0001) %>%
  filter(!is.na(entrezgene_id)) %>%
  filter(entrezgene_id != ensembl_gene_id) %>%
  filter(!grepl(pattern = "MT\\-", x = hgnc_symbol))

# prepare inputs for clusterprofiler -----------------------------------------------------
entrezgene_ids_test <- unique(genes_test_df$entrezgene_id)
entrezgene_ids_universe <- unique(genes_background_df$entrezgene_id)

# test over-representation analysis and gene set enrichment using wikipathway ------------------------------
enricher_out <- tryCatch(expr = enricher(gene = entrezgene_ids_test, TERM2GENE = wp2gene, pvalueCutoff = 1, universe = entrezgene_ids_universe),
                         error = function(e) {warning("ORA failed.");return(NULL)})

if (length(enricher_out) > 0 ) {
  ## convert entrez id to gene symbol
  enricher_out <- setReadable(enricher_out, org.Hs.eg.db, keyType = "ENTREZID")
  enricher_out_all_df <- enricher_out@result
}

# plot enrichment map -----------------------------------------------------
p <- dotplot(object = enricher_out, showCategory=min(50, nrow(enricher_out_all_df[enricher_out_all_df$p.adjust < 0.05,])))
file2write <- paste(dir_out, "dotplot.pdf")
pdf(file2write, width = 12, height = 9, useDingbats = F)
print(p)
dev.off()

enricher_out_pairwise <- enrichplot::pairwise_termsim(enricher_out)
p <- emapplot(x = enricher_out_pairwise,showCategory = min(50, nrow(enricher_out_all_df[enricher_out_all_df$p.adjust < 0.05,]))) 
file2write <- paste(dir_out, "emapplot.pdf")
pdf(file2write, width = 15, height = 10, useDingbats = F)
print(p)
dev.off()

# save output -------------------------------------------------------------
# store results
file2write <- paste0(dir_out, "ORA_Results", ".RDS")
saveRDS(object = enricher_out, file = file2write, compress = T)
## store results
file2write <- paste0(dir_out, "ORA_Results", ".tsv")
write.table(x = enricher_out_all_df, file = file2write, quote = F, row.names = F, sep = "\t")
