# Yige Wu @WashU Jun 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(clusterProfiler)
library(org.Hs.eg.db)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input wikipathway 
wp2gene <- read.gmt("./Resources/Knowledge/Databases/MSigDB/h.all.v7.4.entrez.gmt")
## input degs
daps_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/bap1/annotate_peaks/annotate_BAP1_vs_NonMutant_daps/20210625.v1/BAP1_vs_NonMutant_DAP2Gene.EnhancerPromoter.20210625.v1.tsv")
## input background genes
peaks_background_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/annotate_peaks/annotate_all_peaks_promoters_enhancers/20210625.v1/ccRCC_vs_PT_DAP2Gene.EnhancerPromoter.20210625.v1.tsv") 

# convert gene symbol to entrez ids ---------------------------------------
## add entrez ids to the deg table
daps_df$entrezgene_id <- mapvalues(x = daps_df$Gene, from = peaks_background_df$Gene, to = as.vector(peaks_background_df$entrezgene_id))

genes_process_df <- daps_df %>%
  filter(DAP_direction == "Down") %>%
  filter(peak2gene_type == "Promoter") %>%
  mutate(gene = Gene) %>%
  filter(!is.na(entrezgene_id)) %>%
  filter(entrezgene_id != gene) %>%
  filter(!grepl(pattern = "MT\\-", x = gene))

# prepare inputs for clusterprofiler -----------------------------------------------------
entrezgene_ids_test <- unique(genes_process_df$entrezgene_id)
entrezgene_ids_universe <- unique(peaks_background_df$entrezgene_id[peaks_background_df$peak2gene_type == "Promoter"])

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
pdf(file2write, width = 7, height = 4, useDingbats = F)
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
