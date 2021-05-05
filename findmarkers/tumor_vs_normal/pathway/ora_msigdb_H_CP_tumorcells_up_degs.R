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
version_tmp <- 2
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
deg_all_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/summarize_deg/summarize_tumor_vs_pt_DEGs/20210429.v1/Tumor_DEGs.EnoughDataPoints.Consistent.20210429.v1.tsv")

# convert gene symbol to entrez ids ---------------------------------------
deg_df <- deg_all_df %>%
  filter(Tumor_vs_PT == "Up") %>%
  mutate(gene = genesymbol_deg) %>%
  mutate(avg_logFC = mean_avg_logFC) %>%
  arrange(desc(avg_logFC))
genes2convert <- unique(deg_df$gene)
## retrieve entrezgene_id
genesymbol2entrezid_df <- getBM(attributes=c('entrezgene_id', 'hgnc_symbol'), 
                                filters = 'hgnc_symbol', 
                                values = genes2convert, 
                                mart = ensembl)
## add entrez ids to the deg table
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
# test over-representation analysis and gene set enrichment using wikipathway ------------------------------
enricher_out <- tryCatch(expr = enricher(gene = names(de_genes), TERM2GENE = wp2gene, pvalueCutoff = 1),
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

enricher_out_pairwise <- enrichplot::pairwise_termsim(enricher_out)
p <- emapplot(x = enricher_out_pairwise,showCategory = nrow(enricher_out_all_df[enricher_out_all_df$p.adjust < 0.05,])) 
file2write <- paste(dir_out, "emapplot.pdf")
pdf(file2write, width = 15, height = 10, useDingbats = F)
print(p)
dev.off()

# remove redundancy -------------------------------------------------------
file2write <- paste0(dir_out, "remove_redundancy.txt")
sink(file = file2write)
enricher_out_df <- enricher_out_all_df %>%
  filter(p.adjust < 0.05)
cutoff_frac_overlap <- 0.5
genesetids_test <- enricher_out_df$ID
for (i in genesetids_test[1]) {
  if (i %in% enricher_out_df$ID) {
    genes_i <- str_split(string = enricher_out_df[i, "geneID"],pattern = "\\/")[[1]]
    rowid_i <- which(rownames(enricher_out_df) == i)
    genesetids_j <- enricher_out_df$ID
    for (j in genesetids_j) {
      genes_j <- str_split(string = enricher_out_df[j, "geneID"],pattern = "\\/")[[1]]
      genes_int_ij <- intersect(genes_i, genes_j)
      if (length(genes_int_ij)/length(genes_i) >= cutoff_frac_overlap | length(genes_int_ij)/length(genes_j) >= cutoff_frac_overlap) {
        genes_i_diff_j <- setdiff(x = genes_i, genes_j)
        enricher_out_i_diff_j <- tryCatch(expr = enricher(gene = deg_filtered_df$entrezgene_id[deg_filtered_df$gene %in% genes_i_diff_j], TERM2GENE = wp2gene, pvalueCutoff = 0.05),
                                          error = function(e) {warning("ORA failed.");return(NULL)})
        
        genes_j_diff_i <- setdiff(genes_j, genes_i)
        enricher_out_j_diff_i <- tryCatch(expr = enricher(gene = deg_filtered_df$entrezgene_id[deg_filtered_df$gene %in% genes_j_diff_i], TERM2GENE = wp2gene, pvalueCutoff = 0.05),
                                          error = function(e) {warning("ORA failed.");return(NULL)})
        if (is.null(enricher_out_i_diff_j) & is.null(enricher_out_j_diff_i)) {
          if (length(genes_i) >= length(genes_j)) {
            enricher_out_df <- enricher_out_df[!(rownames(enricher_out_df) %in% j),]
            cat(paste0("<", j, "> Removed because of overlapping with <", i, ">\n\n"))
          } else {
            enricher_out_df <- enricher_out_df[!(rownames(enricher_out_df) %in% i),]
            cat(paste0("<", i, "> Removed because of overlapping with <", j, ">\n\n"))
            next()
          }
        } else if (!is.null(enricher_out_i_diff_j)) {
          cat(paste0("<", i, "> Showed enrichment independent to <", j, ">\n"))
          if (!is.null(enricher_out_j_diff_i)) {
            cat(paste0("<", j, "> Showed enrichment independent to <", i, ">\n"))
            cat("Keep both\n")
            cat("\n")
          } else {
            enricher_out_df <- enricher_out_df[!(rownames(enricher_out_df) %in% j),]
            cat(paste0("<", j, "> Removed because of overlapping with <", i, ">\n\n"))
          }
        } else {
          enricher_out_df <- enricher_out_df[!(rownames(enricher_out_df) %in% i),]
          cat(paste0("<", i, "> Removed because of overlapping with <", j, ">\n\n"))
          next()
        }
      }
    }
  } else {
    next()
  }
}
sink()

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
