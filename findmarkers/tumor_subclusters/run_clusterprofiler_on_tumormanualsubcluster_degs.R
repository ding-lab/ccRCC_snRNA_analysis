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
## input DEGs
markers_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/findallmarkers_wilcox_tumorcells_by_manualsubcluster/20200427.v1/Tumormanualsubcluster.FindAllMarkers.Wilcox.Minpct0.1.Logfc0.25.tsv")
nrow(markers_df)
## filter markers
markers_df <- markers_df %>%
  filter(p_val_adj < 0.05)
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200505.v1/meta_data.20200505.v1.tsv")
## add cluster name 
markers_df$id_aliquot_wu <- mapvalues(x = markers_df$id_aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
markers_df <- markers_df %>%
  mutate(name_cluster = paste0(id_aliquot_wu, "_C", (cluster + 1)))

# convert gene symbol to entrez ids ---------------------------------------
## get biomart
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)
## retrieve UCSC stable ids
genesymbol2ucsc_df <- getBM(attributes=c('entrezgene_id', 'hgnc_symbol'), 
                            filters = 'hgnc_symbol', 
                            values = unique(markers_df$gene), 
                            mart = ensembl)
## add entrez ids to the deg table
markers_df$entrezgene_id <- mapvalues(x = markers_df$gene, from = genesymbol2ucsc_df$hgnc_symbol, to = as.vector(genesymbol2ucsc_df$entrezgene_id))
## filter markers by entrez id
markers_df <- markers_df %>%
  filter(!is.na(entrezgene_id)) %>%
  filter(entrezgene_id != gene)
nrow(markers_df)

# process gmt -------------------------------------------------------------
wp2gene <- read.gmt("./Resources/Databases/Wikipathways/wikipathways-20200510-gmt-Homo_sapiens.gmt")
wp2gene <- wp2gene %>% tidyr::separate(ont, c("name","version","wpid","org"), "%")
wpid2gene <- wp2gene %>% dplyr::select(wpid, gene) #TERM2GENE
wpid2name <- wp2gene %>% dplyr::select(wpid, name) #TERM2NAME

# make ranked gene list ---------------------------------------------------
## initiate result 
ovaresult_df <- NULL
ovaresult_list <- list()
ovaresult_list[["up"]] <- list()
ovaresult_list[["down"]] <- list()
gsearesult_df <- NULL
gsearesult_list <- list()
for (name_cluster_tmp in unique(markers_df$name_cluster)) {
  # name_cluster_tmp <- "VHL"
  ## get genes
  markers_for_cluster_df <- markers_df %>%
    filter(name_cluster == name_cluster_tmp)
  markers_for_cluster_df <- as.data.frame(markers_for_cluster_df)
  ## get gene list
  de_genes <- markers_for_cluster_df$avg_logFC
  de_genes
  names(de_genes) <- markers_for_cluster_df$entrezgene_id
  ## sort
  de_genes <- sort(de_genes, decreasing = T)
  de_genes
  
  gsea_tmp <- tryCatch(expr = GSEA(geneList = de_genes, TERM2GENE = wpid2gene, TERM2NAME = wpid2name, verbose=FALSE, pvalueCutoff = 0.1),
                   error = function(e) {warning("GSEA failed.");return(NULL)})
  if (length(gsea_tmp) > 0) {
    gsea_tmp <- setReadable(gsea_tmp, org.Hs.eg.db, keyType = "ENTREZID")
    gsearesult_list[[name_cluster_tmp]] <- gsea_tmp
    gsea_tmp_df <- as.data.frame(gsea_tmp)
    gsea_tmp_df <- gsea_tmp_df %>%
      mutate(name_cluster = name_cluster_tmp)
    gsearesult_df <- rbind(gsearesult_df, gsea_tmp_df)
  }
  
  # test over-representation analysis and gene set enrichment using wikipathway ------------------------------
  ora_up_tmp <- tryCatch(expr = enricher(gene = names(de_genes[de_genes>0]), TERM2GENE = wpid2gene, TERM2NAME = wpid2name, pvalueCutoff = 0.1),
                   error = function(e) {warning("ORA failed.");return(NULL)})
  if (length(ora_up_tmp) > 0 ) {
    ora_up_tmp <- setReadable(ora_up_tmp, org.Hs.eg.db, keyType = "ENTREZID")
    ovaresult_list[["up"]][[name_cluster_tmp]] <- ora_up_tmp
    ora_up_tmp_df <- as.data.frame(ora_up_tmp)
    ora_up_tmp_df <- ora_up_tmp_df %>%
      mutate(name_cluster = name_cluster_tmp) %>%
      mutate(direction_foldchange = "up")
    ovaresult_df <- rbind(ovaresult_df, ora_up_tmp_df)
  }
  
  # test over-representation analysis and gene set enrichment using wikipathway ------------------------------
  ora_down_tmp <- tryCatch(expr = enricher(gene = names(de_genes[de_genes<0]), TERM2GENE = wpid2gene, TERM2NAME = wpid2name, pvalueCutoff = 0.1),
                           error = function(e) {warning("ORA failed.");return(NULL)})
  if (length(ora_down_tmp) > 0) {
    ora_down_tmp <- setReadable(ora_down_tmp, org.Hs.eg.db, keyType = "ENTREZID")
    ovaresult_list[["up"]][[name_cluster_tmp]] <- ora_down_tmp
    ora_down_tmp_df <- as.data.frame(ora_down_tmp)
    ora_down_tmp_df <- ora_down_tmp_df %>%
      mutate(name_cluster = name_cluster_tmp) %>%
      mutate(direction_foldchange = "down")
    ovaresult_df <- rbind(ovaresult_df, ora_down_tmp_df)
  }
}

# save outputs ------------------------------------------------------------
## save data frames
file2write <- paste0(dir_out, "ORA.", "Result.", run_id, ".tsv")
write.table(x = ovaresult_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "GSEA.", "Result.", run_id, ".tsv")
write.table(x = gsearesult_df, file = file2write, quote = F, sep = "\t", row.names = F)
## save lists
file2write <- paste0(dir_out, "ORA.", "Result.", run_id, ".RDS")
saveRDS(file = file2write, object = ovaresult_list, compress = T)
file2write <- paste0(dir_out, "GSEA.", "Result.", run_id, ".RDS")
saveRDS(file = file2write, object = gsearesult_list, compress = T)


