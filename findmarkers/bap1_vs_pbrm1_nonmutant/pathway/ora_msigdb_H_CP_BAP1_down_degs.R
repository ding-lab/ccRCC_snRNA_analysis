# Yige Wu @WashU Sep 2020

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
wp2gene1 <- read.gmt("./Resources/Knowledge/Databases/MSigDB/h.all.v7.4.entrez.gmt")
wp2gene2 <- read.gmt("./Resources/Knowledge/Databases/MSigDB/c2.cp.v7.4.entrez.gmt")
wp2gene <- rbind(wp2gene1, wp2gene2)
## input degs
# deg_all_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/bap1_vs_pbrm1_nonmutant/summarize_degs/summarize_BAP1_vs_PBRM1_NonMutant_DEGs/20210430.v1/BAP1_DEGs.Consistent20210430.v1.tsv")
deg_all_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/bap1_vs_pbrm1_nonmutant/summarize_degs/unite_BAP1_snRNA_bulkRNA_protein_DEGs/20210610.v1/BAP1_snRNA_DEGs.Consistent.CNVcorrected.20210610.v1.tsv")
genes_background_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/bap1_vs_pbrm1_nonmutant/findallmarker_LR_all_BAP1_tumorcells_vs_PBRM1_NonMutant_cells_on_katmai/20210609.v1/gProfiler_hsapiens_6-10-2021_5-56-28 PM.csv")

# convert gene symbol to entrez ids ---------------------------------------
genes2convert <- unique(genes_background_df$genesymbol_deg)
genes2convert <- unique(deg_df$gene)

## retrieve entrezgene_id
genesymbol2entrezid_df <- getBM(attributes=c('entrezgene_id', 'hgnc_symbol'), 
                                filters = 'hgnc_symbol', 
                                values = genes2convert, 
                                mart = ensembl)
## add entrez ids to the deg table
deg_df <- deg_all_df %>%
  filter(BAP1_vs_OtherTumor_snRNA == "Down") %>%
  mutate(gene = genesymbol_deg) %>%
  mutate(avg_logFC = mean_avg_logFC) %>%
  arrange(desc(avg_logFC))
deg_df$entrezgene_id <- mapvalues(x = deg_df$gene, from = genesymbol2entrezid_df$hgnc_symbol, to = as.vector(genesymbol2entrezid_df$entrezgene_id))
## filter markers by entrez id
deg_filtered_df <- deg_df %>%
  filter(!is.na(entrezgene_id)) %>%
  filter(entrezgene_id != gene) %>%
  filter(!grepl(pattern = "MT\\-", x = gene))

# prepare inputs for clusterprofiler -----------------------------------------------------
de_genes <- deg_filtered_df$avg_logFC
de_genes
names(de_genes) <- deg_filtered_df$entrezgene_id
## sort
de_genes <- sort(de_genes, decreasing = T)
de_genes
## make universe gene list
universe_genes <- as.character(genes_background_df$converted_alias[!is.na(genes_background_df$converted_alias)])
# test over-representation analysis and gene set enrichment using wikipathway ------------------------------
enricher_out <- tryCatch(expr = enricher(gene = names(de_genes), TERM2GENE = wp2gene, pvalueCutoff = 1, universe = universe_genes),
                         error = function(e) {warning("ORA failed.");return(NULL)})

if (length(enricher_out) > 0 ) {
  ## convert entrez id to gene symbol
  enricher_out <- setReadable(enricher_out, org.Hs.eg.db, keyType = "ENTREZID")
  enricher_out_all_df <- enricher_out@result
}

# plot enrichment map -----------------------------------------------------
p <- dotplot(object = enricher_out, showCategory=min(50, nrow(enricher_out_all_df[enricher_out_all_df$p.adjust < 0.05,])))
file2write <- paste(dir_out, "dotplot.pdf")
pdf(file2write, width = 12, height = 10, useDingbats = F)
print(p)
dev.off()

enricher_out_pairwise <- enrichplot::pairwise_termsim(enricher_out)
p <- emapplot(x = enricher_out_pairwise,showCategory = min(50, nrow(enricher_out_all_df[enricher_out_all_df$p.adjust < 0.05,]))) 
file2write <- paste(dir_out, "emapplot.pdf")
pdf(file2write, width = 20, height = 15, useDingbats = F)
print(p)
dev.off()

# save output -------------------------------------------------------------
## store results
file2write <- paste0(dir_out, "Fold_Changes", ".RDS")
saveRDS(object = de_genes, file = file2write, compress = T)
## store results
file2write <- paste0(dir_out, "ORA_Results", ".RDS")
saveRDS(object = enricher_out, file = file2write, compress = T)
## store results
file2write <- paste0(dir_out, "ORA_Results", ".tsv")
write.table(x = enricher_out_all_df, file = file2write, quote = F, row.names = F, sep = "\t")
