# Yige Wu @WashU Sep 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input summary for cell-cell interactin
summary_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/summarize_interactions/rank_across_celltypes_by_pair_by_avgsigmean/20200925.v3/cellphonedb.summary_across_celltypes_by_pair.min5.tsv")
## input the average expression calculated (SCT)
avgexp_df <- fread(input = "./Resources/Analysis_Results/average_expression/format_expression/format_avgexp_sct_usedata_bycelltypeshorter/20200907.v2/formated.averageexpression.SCT.slotdata.bycelltype.shorter.31_aliquot_integration.20200907.v2.tsv", data.table = F)
## input complex annotation
complex2gene_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/annotate_interactions/annotate_cellphonedb_complex_to_genesymbols/20200928.v1/CellPhoneDB_complexs_genesymbols.tsv")

# get all gene names from cellphone db ------------------------------------
genes_cp_all <- unique(c(unique(summary_df$gene.source), unique(summary_df$gene.target)))
length(genes_cp_all)
## convert the complex to gene names
genes_complexesmapped <- mapvalues(x = genes_cp_all, from = complex2gene_df$identifiers_CellPhoneDB, to = as.vector(complex2gene_df$components_genesymbols))
genes_complexesmapped
## split
genes_complexessplit_list <- sapply(X = genes_complexesmapped, FUN = function(x) {
  text_split <- str_split(string = x, pattern = "_")[[1]]
  return(text_split)
})
genes_complexessplit_list
names(genes_complexessplit_list) <- genes_cp_all
genes_complexessplit_vec <- unique(unlist(genes_complexessplit_list))
length(genes_complexessplit_vec)

# filter average expression and get the top cell type ---------------------
avgexp_filtered_df <- avgexp_df %>%
  filter(gene %in% genes_complexessplit_vec)
avgexp_melt_df <- melt(data = avgexp_filtered_df)
avgexp_ordered_df <- avgexp_melt_df %>%
  group_by(gene) %>%
  dplyr::mutate(rank_singlegene_acrosscelltypes = order(order(value, decreasing = T))) %>%
  rename(Celltype.colname = variable)

# make cell type name mapping ---------------------------------------------
celltype2colname_df <- data.frame(Celltype = unique(c(summary_df$Cell_type.source, summary_df$Cell_type.target)))
celltype2colname_df <- celltype2colname_df %>%
  mutate(Celltype.colname = gsub(x = Celltype, pattern = ' |\\/|\\-', replacement = "."))
## add cell type colname
avgexp_ordered_df$Celltype <- mapvalues(x = avgexp_ordered_df$Celltype.colname, from = celltype2colname_df$Celltype.colname, to = as.vector(celltype2colname_df$Celltype))

# map the rank of single ligand/receptor gene exp -------------------------
summary_df$rank_genesource_acrosscelltypes <- sapply(X = 1:nrow(summary_df), FUN = function(x, complex2gene_list, rank_df, data_df) {
  gene_tmp = data_df[x, "gene.source"]
  celltype_tmp = data_df[x, "Cell_type.source"]
  genesymbols <- complex2gene_list[[gene_tmp]]
  ranks_tmp <- rank_df$rank_singlegene_acrosscelltypes[rank_df$gene %in% genesymbols & rank_df$Celltype == celltype_tmp]
  rank_min = min(ranks_tmp)
  return(rank_min)
}, complex2gene_list = genes_complexessplit_list, rank_df = avgexp_ordered_df, data_df = summary_df)
summary_df$rank_genetarget_acrosscelltypes <- sapply(X = 1:nrow(summary_df), FUN = function(x, complex2gene_list, rank_df, data_df) {
  gene_tmp = data_df[x, "gene.target"]
  celltype_tmp = data_df[x, "Cell_type.target"]
  genesymbols <- complex2gene_list[[gene_tmp]]
  ranks_tmp <- rank_df$rank_singlegene_acrosscelltypes[rank_df$gene %in% genesymbols & rank_df$Celltype == celltype_tmp]
  rank_min = min(ranks_tmp)
  return(rank_min)
}, complex2gene_list = genes_complexessplit_list, rank_df = avgexp_ordered_df, data_df = summary_df)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "cellphonedb.summary_across_celltypes_by_pair.min5.singlegeneranked.tsv")
write.table(x = summary_df, file = file2write, sep = "\t", quote = F, row.names = F)
file2write <- paste0(dir_out, "cellphonedb.complex2gene.RDS")
saveRDS(object = genes_complexessplit_list, file = file2write, compress = T)

