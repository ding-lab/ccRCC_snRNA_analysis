# Yige Wu @WashU Apr 2020

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

# input denpendencies -----------------------------------------------------
## input average expression by cell type by aliquot
avgexp_bycelltype_byaliquot_df <- fread(input = "./Resources/Analysis_Results/average_expression/averageexpression_bycelltypeshorter_byaliquot_on_katmai/20200411.v1/averageexpression_bycelltypeshorter.30_aliquot_integration.20200411.v1.tsv", data.table = F)
## input average expression by aliquot
avgexp_byaliquot_df <- fread(input = "./Resources/Analysis_Results/average_expression/averageexpression_byaliquot_on_katmai/20200411.v1/averageexpression_byaliquot.30_aliquot_integration.20200411.v1.tsv", data.table = F)

# make adjustment for the average expression by celltype and by aliquot---------------------------------------------------------
## get unique filtered genes
genes_filtered <- unique(celltypemarker_filtered_df$gene)
length(genes_filtered)
## filter average expression by filtered genes
avgexp_bycelltype_byaliquot_filtered_df <- avgexp_bycelltype_byaliquot_df %>%
  filter(V1 %in% genes_filtered)
## melt
avgexp_bycelltype_byaliquot_melt_df <- melt(data = avgexp_bycelltype_byaliquot_filtered_df)
avgexp_bycelltype_byaliquot_melt_df <- avgexp_bycelltype_byaliquot_melt_df %>%
  rename(gene_symbol = V1) %>%
  mutate(id_aliquot_celltype = str_split_fixed(string = variable, pattern = "\\.", n = 2)[,2]) %>%
  mutate(aliquot = str_split_fixed(string = id_aliquot_celltype, pattern = "_", n = 2)[,1]) %>%
  mutate(celltype = str_split_fixed(string = id_aliquot_celltype, pattern = "_", n = 2)[,2]) %>%
  rename(avg_exp = value) %>%
  select(gene_symbol,aliquot, celltype, avg_exp)

# make adjustment for the average expression by celltype and by aliquot---------------------------------------------------------
## filter average expression by filtered genes
avgexp_byaliquot_filtered_df <- avgexp_byaliquot_df %>%
  filter(V1 %in% genes_filtered)
## melt
avgexp_byaliquot_melt_df <- melt(data = avgexp_byaliquot_filtered_df)
avgexp_byaliquot_melt_df <- avgexp_byaliquot_melt_df %>%
  rename(gene_symbol = V1) %>%
  mutate(aliquot = str_split_fixed(string = variable, pattern = "\\.", n = 2)[,2]) %>%
  mutate(celltype = "All_Cells") %>%
  rename(avg_exp = value) %>%
  select(gene_symbol, aliquot, celltype, avg_exp)

# merge -------------------------------------------------------------------
avgexp_df <- rbind(avgexp_bycelltype_byaliquot_melt_df, avgexp_byaliquot_melt_df)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "avgexp_bycelltype.", "long_data_frame.", run_id, ".tsv")
write.table(x = avgexp_df, file = file2write, quote = F, sep = "\t", row.names = F)

