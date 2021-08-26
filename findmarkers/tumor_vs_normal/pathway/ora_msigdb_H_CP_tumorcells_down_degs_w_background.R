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
# deg_all_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/summarize_deg/summarize_tumor_vs_pt_DEGs/20210419.v1/Tumor_DEGs.EnoughDataPoints.Consistent.20210419.v1.tsv")
# deg_all_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/summarize_deg/summarize_tumor_vs_pt_DEGs/20210429.v1/Tumor_DEGs.EnoughDataPoints.Consistent.20210429.v1.tsv")
deg_all_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/summarize_deg/unite_tumor_vs_normal_snRNA_individual_and_CNVcorrected_DEGs/20210824.v1/Consistent.Tumor_vs_PT_DEGs.CNVcorrected.20210824.v1.tsv")
genes_background_df <- fread(data.table = F,input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/findmarker_LR_all_ccRCC_vs_pt_on_katmai/20210824.v1/LR.logfc.threshold0.min.pct0.1.min.diff.pct0.AssayRNA.tsv")

# convert gene symbol to entrez ids ---------------------------------------
genes2convert <- unique(genes_background_df$genesymbol_deg)
## retrieve entrezgene_id
genesymbol2entrezid_df <- getBM(attributes=c('entrezgene_id', 'hgnc_symbol'), 
                                filters = 'hgnc_symbol', 
                                values = genes2convert, 
                                mart = ensembl)
## filter DEGs
deg_df <- deg_all_df %>%
  filter(!is.na(FDR.CNVcorrected) & FDR.CNVcorrected < 0.05) %>%
  filter(Num_sig_up == 0) %>%
  mutate(gene = genesymbol_deg) %>%
  mutate(avg_log2FC = avg_log2FC.allTumorcellsvsPT) %>%
  arrange(desc(avg_log2FC))

## add entrez ids to the deg table
deg_df$entrezgene_id <- mapvalues(x = deg_df$gene, from = genesymbol2entrezid_df$hgnc_symbol, to = as.vector(genesymbol2entrezid_df$entrezgene_id))
## filter markers by entrez id
deg_filtered_df <- deg_df %>%
  filter(!is.na(entrezgene_id)) %>%
  filter(entrezgene_id != gene) %>%
  filter(!grepl(pattern = "MT\\-", x = gene))

# prepare inputs for clusterprofiler -----------------------------------------------------
de_genes <- deg_filtered_df$avg_log2FC
de_genes
names(de_genes) <- deg_filtered_df$entrezgene_id
## sort
de_genes <- sort(de_genes, decreasing = T)
de_genes
# test over-representation analysis and gene set enrichment using wikipathway ------------------------------
background_genes <- genesymbol2entrezid_df$entrezgene_id[!is.na(genesymbol2entrezid_df$entrezgene_id)]; background_genes <- as.character(background_genes)
enricher_out <- tryCatch(expr = enricher(gene = names(de_genes), TERM2GENE = wp2gene, pvalueCutoff = 1, universe = background_genes),
                         error = function(e) {warning("ORA failed.");return(NULL)})

if (length(enricher_out) > 0 ) {
  ## convert entrez id to gene symbol
  enricher_out <- setReadable(enricher_out, org.Hs.eg.db, keyType = "ENTREZID")
  enricher_out_all_df <- enricher_out@result
}

# plot enrichment map -----------------------------------------------------
p <- dotplot(object = enricher_out, showCategory=nrow(enricher_out_all_df[enricher_out_all_df$p.adjust < 0.05,]))
file2write <- paste(dir_out, "dotplot.pdf")
pdf(file2write, width = 12, height = 6, useDingbats = F)
print(p)
dev.off()

enricher_out_pairwise <- enrichplot::pairwise_termsim(x = enricher_out, showCategory = 50)
View(enricher_out_pairwise@termsim)
p <- emapplot(x = enricher_out_pairwise,showCategory = nrow(enricher_out_all_df[enricher_out_all_df$p.adjust < 0.05,])) 
file2write <- paste(dir_out, "emapplot.pdf")
pdf(file2write, width = 15, height = 10, useDingbats = F)
print(p)
dev.off()

# save output -------------------------------------------------------------
## store results
file2write <- paste0(dir_out, "Fold_Changes", ".RDS")
saveRDS(object = de_genes, file = file2write, compress = T)
## store results
file2write <- paste0(dir_out, "ORA_Results", ".RDS")
saveRDS(object = enricher_out_pairwise, file = file2write, compress = T)
## store results
file2write <- paste0(dir_out, "ORA_Results", ".tsv")
write.table(x = enricher_out_all_df, file = file2write, quote = F, row.names = F, sep = "\t")
