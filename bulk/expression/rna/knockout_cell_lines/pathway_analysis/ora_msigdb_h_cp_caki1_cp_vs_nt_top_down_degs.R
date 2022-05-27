# Yige Wu @WashU May 2022

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
setwd(dir_base)
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "clusterProfiler",
  "org.Hs.eg.db",
  "ggnewscale"
)
for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    print(paste0(pkg_name_tmp, "is being installed!"))
    BiocManager::install(pkgs = pkg_name_tmp, update = F)
    install.packages(pkg_name_tmp, dependencies = T)
  }
  print(paste0(pkg_name_tmp, " is installed!"))
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input -------------------------------------------------------------------
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/expression/rna/knockout_cell_lines/DEG_analysis/run_edgeR_DEG_caki_cp_vs_caki_NT/20220517.v2/Caki1_CP_vs_Caki1_NT.DEGs.20220517.v2.tsv")
## input gene-set data
gs2gene1 <- read.gmt("./Resources/Knowledge/Databases/MSigDB/msigdb_v7.4_GMTs/h.all.v7.4.entrez.gmt")
gs2gene2 <- read.gmt("./Resources/Knowledge/Databases/MSigDB/msigdb_v7.4_GMTs/c2.cp.v7.4.entrez.gmt")
gs2gene <- rbind(gs2gene1, gs2gene2)

# prepare data for analysis -----------------------------------------------
genes_background_df <- deg_df %>%
  filter(!is.na(entrezgene) & gene_biotype == "protein_coding")
nrow(genes_background_df)
genes_test_df <- genes_background_df %>%
  filter(logFC < -1 & FDR < 0.05) %>%
  arrange(logFC)
# entrezgene_ids_test <- head(unique(genes_test_df$entrezgene), 500);
entrezgene_ids_test <- head(unique(genes_test_df$entrezgene), 1000);
entrezgene_ids_test <- as.character(entrezgene_ids_test); length(entrezgene_ids_test)
entrezgene_ids_background <- unique(genes_background_df$entrezgene); entrezgene_ids_background <- as.character(entrezgene_ids_background); length(entrezgene_ids_background)
genes_test_df %>%
  filter(entrezgene %in% entrezgene_ids_test) %>%
  View()

# test over-representation analysis and gene set enrichment ------------------------------
enricher_out <- tryCatch(expr = enricher(gene = entrezgene_ids_test, TERM2GENE = gs2gene, pvalueCutoff = 1, universe = entrezgene_ids_background),
                         error = function(e) {warning("ORA failed.");return(NULL)})
if (length(enricher_out) > 0 ) {
  ## convert entrez id to gene symbol
  enricher_out <- setReadable(enricher_out, org.Hs.eg.db, keyType = "ENTREZID")
  enricher_out_all_df <- enricher_out@result
}

# plot enrichment map -----------------------------------------------------
enricher_out_pairwise <- enrichplot::pairwise_termsim(enricher_out)
p <- emapplot(x = enricher_out_pairwise,showCategory = min(50, nrow(enricher_out_all_df[enricher_out_all_df$p.adjust < 0.05,]))) 
file2write <- paste(dir_out, "emapplot.pdf")
pdf(file2write, width = 25, height = 20, useDingbats = F)
print(p)
dev.off()

# save output -------------------------------------------------------------
## store results
file2write <- paste0(dir_out, "ORA_Results", ".tsv")
write.table(x = enricher_out_all_df, file = file2write, quote = F, row.names = F, sep = "\t")

