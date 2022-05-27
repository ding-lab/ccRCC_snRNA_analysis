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
## input gene-set data
gs2gene1 <- read.gmt("./Resources/Knowledge/Databases/MSigDB/h.all.v7.4.entrez.gmt")
gs2gene2 <- read.gmt("./Resources/Knowledge/Databases/MSigDB/c2.cp.v7.4.entrez.gmt")
gs2gene <- rbind(gs2gene1, gs2gene2)
## input degs
dep_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/protein/compare_bulk_protein_BAP1_tumor_vs_NonMutant_tumors/20210913.v1/Bulk_Protein_BAP1_vs_NonMutant.Wilcox.20210913.v1.tsv")

# convert gene symbol to entrez ids ---------------------------------------
genes_background_df <- dep_df %>%
  filter(number_group1 >= 10 & number_group2 >= 10)
genes2convert <- unique(genes_background_df$gene_symbol)
## retrieve entrezgene_id
ensembl <- useMart("ensembl")
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
genesymbol2entrezid_df <- getBM(attributes=c('entrezgene_id', 'hgnc_symbol'), 
                                filters = 'hgnc_symbol', 
                                values = genes2convert, 
                                mart = ensembl)
## add entrez ids to the deg table
genes_background_df$entrezgene_id <- mapvalues(x = genes_background_df$gene, from = genesymbol2entrezid_df$hgnc_symbol, to = as.vector(genesymbol2entrezid_df$entrezgene_id))
##
genes_test_df <- genes_background_df %>%
  filter(meddiff_exp < 0) %>%
  filter(!is.na(entrezgene_id)) %>%
  filter(entrezgene_id != gene_symbol) %>%
  filter(!grepl(pattern = "MT\\-", x = gene_symbol))

# prepare inputs for clusterprofiler -----------------------------------------------------
entrezgene_ids_test <- unique(genes_test_df$entrezgene_id)
entrezgene_ids_background <- unique(genes_background_df$entrezgene_id)

# test over-representation analysis and gene set enrichment using wikipathway ------------------------------
enricher_out <- tryCatch(expr = enricher(gene = entrezgene_ids_test, TERM2GENE = gs2gene, pvalueCutoff = 1, universe = entrezgene_ids_background),
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
