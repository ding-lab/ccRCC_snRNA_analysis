# Yige Wu @WashU Jun 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library(biomaRt)
library(org.Hs.eg.db)
library(dplyr)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input peak fold changes
peaks_anno_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Peak_Annotation/28_snATACmerged_allPeaks.Annotated.20210712.tsv")
## input coaccessiblity results
coaccess_peak2genes_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/annotate_peaks/annotate_coaccessible_peaks_28samples/20211011.v1/Coaccessible_Peaks.Annotated.20211011.v1.tsv")
## get biomart
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)
filters = listFilters(ensembl)

# annotate and filter peaks ------------------------------------------------------------
## annotate peaks
peaks_anno_df <- peaks_anno_df %>%
  dplyr::mutate(peak2gene_type = str_split_fixed(string = Type, pattern = " \\(", n = 2)[,1]) %>%
  # rename(Gene = SYMBOL) %>%
  dplyr::select(peak, peak2gene_type, Gene)
# peaks_anno_df %>%
#   filter(peak2gene_type == "Promoter") %>%
#   nrow()
# 
# peaks_anno_df %>%
#   filter(peak2gene_type != "Promoter") %>%
#   nrow()

# annotate with coaccessiblity results ------------------------------------
peaks_anno_wcoaccess_df <- merge(x = peaks_anno_df, 
                       y = coaccess_peak2genes_df %>%
                         dplyr::select(Peak1, Peak2, peak2gene_type.2, genesymbol.2, coaccess) %>%
                         dplyr::rename(peak.coaccess = Peak2) %>%
                         dplyr::rename(peak2gene_type.coaccess = peak2gene_type.2) %>%
                         dplyr::rename(genesymbol.coaccess = genesymbol.2) %>%
                         dplyr::rename(coaccess_score = coaccess),
                       by.x = c("peak"), by.y = c("Peak1"))
peak2gene_enhancers_df <- peaks_anno_wcoaccess_df %>%
  dplyr::filter(peak2gene_type != "Promoter" & peak2gene_type.coaccess == "Promoter")
# peak2gene_enhancers_df %>%
#   dplyr::select(peak) %>%
#   unique() %>%
#   nrow()
peak2gene_enh_pro_df <- rbind(peak2gene_enhancers_df %>%
                                mutate(peak2gene_type = "Enhancer") %>%
                                mutate(Gene = genesymbol.coaccess),
                              peaks_anno_df %>%
                                filter(peak2gene_type == "Promoter") %>%
                                mutate(peak.coaccess = NA) %>%
                                mutate(peak2gene_type.coaccess = NA) %>%
                                mutate(genesymbol.coaccess = NA) %>%
                                mutate(coaccess_score = NA))


# annotate with entrez ids ------------------------------------------------
genes2convert <- unique(c(peaks_anno_df$Gene, peak2gene_enh_pro_df$Gene))
## retrieve entrezgene_id
genesymbol2entrezid_df <- getBM(attributes=c('entrezgene_id', 'hgnc_symbol'), 
                                filters = 'hgnc_symbol', 
                                values = genes2convert, 
                                mart = ensembl)
## add entrez ids to the deg table
peaks_anno_df$entrezgene_id <- mapvalues(x = peaks_anno_df$Gene, from = genesymbol2entrezid_df$hgnc_symbol, to = as.vector(genesymbol2entrezid_df$entrezgene_id))
peak2gene_enh_pro_df$entrezgene_id <- mapvalues(x = peak2gene_enh_pro_df$Gene, from = genesymbol2entrezid_df$hgnc_symbol, to = as.vector(genesymbol2entrezid_df$entrezgene_id))

# write outupt ------------------------------------------------------------
file2write <- paste0(dir_out, "ccRCC_vs_PT_DAPs.Annotated.", run_id, ".tsv")
write.table(file = file2write, x = peaks_anno_df, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "ccRCC_vs_PT_DAP2Gene.EnhancerPromoter.", run_id, ".tsv")
write.table(file = file2write, x = peak2gene_enh_pro_df, quote = F, sep = "\t", row.names = F)


