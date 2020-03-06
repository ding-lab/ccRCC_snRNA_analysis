# Yige Wu @WashU Feb 2020
## getting the expressin of CNv-related DEGs in cells with or without CNVs per sample

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## set chr_region to process
chr_region_tmp <- "3p26"
## input case CNV profile
bulk_bicseq_cnv_state_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/bulk/write_sample_bicseq_cnv_profile/20200227.v1/Bulk_WGS_Chr_CNV_Profile.20200227.v1.tsv", data.table = F)
bulk_bicseq_cnv_state_df$`5q`[bulk_bicseq_cnv_state_df$Case == "C3N-00242"] <- "gain"
## input meta data
meta_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/make_meta_data/20191105.v1/meta_data.20191105.v1.tsv", data.table = F)
## input DEG united list
deg_by_chr_region_df <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/copy_number/findallmarkers/unite_degs_by_chr_region/20200228.v1/FindMarkers.Wilcox.ExpectedCNV_vs_Neutral.20200228.v1.tsv", data.table = F)
## input the genetic downstream table
genetic_alt_downstream_genes <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/dependencies/write_ccrcc_genetic_event_downstream_genes/20200302.v2/ccRCC_Genetic_Event_Downstream_Genes.20200302.v2.tsv", data.table = F)

# count reoccurance -------------------------------------------------------
## filter by chromosomal region
deg_by_chr_region_df <- deg_by_chr_region_df %>%
  filter(chr_region == chr_region_tmp)

deg_by_aliquot_count <- deg_by_chr_region_df %>%
  select(gene, aliquot) %>%
  table () %>%
  as.data.frame() %>%
  filter(Freq > 0) %>%
  arrange(desc(Freq))

# filter by case and associated CNV ------------------------------------------------------------
## filter cases by bulk WGS CNV profile
cases_expected_cnv <- bulk_bicseq_cnv_state_df[bulk_bicseq_cnv_state_df$`14q` == "loss",]

## map cases
deg_by_chr_region_df$Case <- mapvalues(x = deg_by_chr_region_df$aliquot, from = meta_tab$Aliquot.snRNA, to = as.vector(meta_tab$Case))

## filter by case
deg_by_chr_region_df <- deg_by_chr_region_df %>%
  filter(Case %in% cases_expected_cnv$Case) %>%
  arrange(p_val_adj, p_val)

# write filtered table for 5q associated DEGs -----------------------------
write.table(x = deg_by_chr_region_df, file = paste0(dir_out, chr_region_tmp, ".FindMarkers.Wilcox.ExpectedCNV_vs_Neutral", ".", run_id, ".tsv"), sep = "\t", quote = F, row.names = F)

# annotate gene to pathways -----------------------------------------------
## inpu over-representation gene set result
ora_result <- fread(input = paste0("Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/copy_number/findallmarkers/examine_", chr_region_tmp, "_degs/20200302.v1/ORA_results.tab"), data.table = F)
## convert ORA table to list, by pathway
path2genes_list <- lapply(X = ora_result$pathway, FUN = function(p, ora_result_df) {
  gene_string <- ora_result_df$members_input_overlap[ora_result_df$pathway == p]
  genes <- str_split(string = gene_string, pattern = "; ")[[1]]
  return(genes)
}, ora_result_df = ora_result)
names(path2genes_list) <- ora_result$pathway
gene2path_df <- data.frame(gene = unlist(path2genes_list), pathway = rep(x = names(path2genes_list), sapply(path2genes_list, FUN = function(gene_list) return(length(gene_list)))))
## write the annotation table
write.table(x = gene2path_df, file = paste0(dir_out, chr_region_tmp, ".FindMarkers.Wilcox.ExpectedCNV_vs_Neutral",  ".", "Gene2Pathway",".", run_id, ".tsv"), sep = "\t", quote = F, row.names = F)

## filter by gene
deg_filtered_by_gene <- deg_by_chr_region_df %>%
  filter(gene %in% genetic_alt_downstream_genes$target_genesymbol)

deg_filtered_by_gene <- merge(deg_filtered_by_gene, genetic_alt_downstream_genes, by.x = c("gene"), by.y = c("target_genesymbol"), all.x = T)



