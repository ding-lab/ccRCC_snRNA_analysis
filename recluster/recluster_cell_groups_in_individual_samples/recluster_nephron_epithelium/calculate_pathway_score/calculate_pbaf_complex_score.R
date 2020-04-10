# Yige Wu @WashU Apr 2020
## running on local
## for calculating VHL-HIF activity score using scaled average expression
## PBAF complex

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
## input scaled average expression
scaled_avgexp_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/other/scale_averageexpression/20200407.v1/averageexpression_tumor_cells_by_manual_subcluster.scaled.20200407.v1.tsv", data.table = F)
## input ccRCC downstream genes
ccrcc_gene_alt_downstream_genes_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/dependencies/write_ccrcc_genetic_event_downstream_genes/20200302.v2/ccRCC_Genetic_Event_Downstream_Genes.20200302.v2.tsv", data.table = F)
## input the variable gene list
variable_genes_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/findvariablefeatures/findvariablefeatures_tumor_cells/20200310.v1/findvariablefeatures_tumor_cells.20200310.v1.tsv", data.table = F)

# filter out rows with all NAs in scaled data--------------------------------------------
scaled_avgexp_df <- scaled_avgexp_df[rowSums(!is.na(scaled_avgexp_df)) > 1,]

# get the list of eligible HIF pathway members ----------------------------
pathway_genes <- pbaf_genes
## filter by genes in scaled data
pathway_genes <- intersect(pathway_genes, as.vector(scaled_avgexp_df$V1))
# pathway_genes <- intersect(pathway_genes, as.vector(variable_genes_df$gene))
pathway_genes

# average by column/manual subcluster -------------------------------------------------------
## filter row by the pathway members
pathway_scaled_avgexp_df <- scaled_avgexp_df %>%
  filter(V1 %in% pathway_genes)
## get matrix-like scaled data
pathway_scaled_avgexp_mat <- pathway_scaled_avgexp_df %>%
  select(-V1)
## average by column
avg_pathway_scaled_avgexp_mat <- colMeans(pathway_scaled_avgexp_mat)
## make long data frame
avg_pathway_scaled_avgexp_df <- data.frame(exp_manual_subcluster_name = str_split_fixed(string = names(avg_pathway_scaled_avgexp_mat), pattern = "\\.", n = 2)[,2],
                                           pathway_score = avg_pathway_scaled_avgexp_mat,
                                           pathway_name = "PBAF")
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "pbaf_complex.", "avg_pathway_scaled_avgexp.", run_id, ".tsv")
write.table(x = avg_pathway_scaled_avgexp_df, file = file2write, quote = F, sep = "\t", row.names = F)

