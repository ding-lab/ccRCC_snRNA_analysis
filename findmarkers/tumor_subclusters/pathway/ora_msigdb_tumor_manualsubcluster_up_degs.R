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
wp2gene <- read.gmt("./Resources/Knowledge/Databases/MSigDB/h.all.v7.4.entrez.gmt")
## input degs
deg_all_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/unite_degs/unite_degs_for_tumor_manualcluster/20210413.v1/TumorManualCluster.DEGs.Wilcox.Minpct0.1.Logfc0.min.diff.pct0.1.tsv")

# convert gene symbol to entrez ids ---------------------------------------
deg_df <- deg_all_df %>%
  filter(avg_logFC > 0 & p_val_adj < 0.001) %>%
  filter(cellcount_group1 >= 50) %>%
  filter(cellcount_group2 >= 50)
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

# prepare inputs for clusterprofiler -----------------------------------------------------
foldchange_list <- list()
ora_result_list <- list()
ora_result_df <- NULL
for (easyid_tmp in unique(deg_df$easy_id)) {
  for (cluster_tmp in unique(deg_df$cluster[deg_df$easy_id == easyid_tmp])) {
    cluster_name_tmp <- paste0(easyid_tmp, "_C", cluster_tmp)
    ## get gene list
    deg_tmp_df <- deg_df %>%
      filter(easy_id == easyid_tmp) %>%
      filter(cluster == cluster_tmp)
    de_genes <- deg_tmp_df$avg_logFC
    de_genes
    names(de_genes) <- deg_tmp_df$entrezgene_id
    ## sort
    de_genes <- sort(de_genes, decreasing = T)
    de_genes
    foldchange_list[[cluster_name_tmp]] <- de_genes
    # test over-representation analysis and gene set enrichment using wikipathway ------------------------------
    enricher_out <- tryCatch(expr = enricher(gene = names(de_genes), TERM2GENE = wp2gene, pvalueCutoff = 1),
                             error = function(e) {warning("ORA failed.");return(NULL)})
    
    if (length(enricher_out) > 0 ) {
      ## convert entrez id to gene symbol
      enricher_out <- setReadable(enricher_out, org.Hs.eg.db, keyType = "ENTREZID")
      ora_result_list[[cluster_name_tmp]] <- enricher_out
      
      enricher_out_df <- enricher_out@result
      ora_result_df <- rbind(ora_result_df, enricher_out_df %>%
                               mutate(easy_id = easyid_tmp) %>%
                               mutate(cluster = cluster_tmp))
    }
    enricher_out <- NULL
  }
}
# save output -------------------------------------------------------------
## store results
file2write <- paste0(dir_out, "Fold_Changes", ".RDS")
saveRDS(object = foldchange_list, file = file2write, compress = T)
## store results
file2write <- paste0(dir_out, "ORA_Results", ".RDS")
saveRDS(object = ora_result_list, file = file2write, compress = T)
## store results
file2write <- paste0(dir_out, "ORA_Results", ".tsv")
write.table(x = ora_result_df, file = file2write, quote = F, row.names = F, sep = "\t")
