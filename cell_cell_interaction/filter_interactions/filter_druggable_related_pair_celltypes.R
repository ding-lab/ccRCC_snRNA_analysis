# Yige Wu @WashU Sep 2020
## reference of the cellphonedb output: https://www.cellphonedb.org/documentation

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
## input the strongest combinations of ligand-receptor pairs and cell types
cellphone_sum_by_pair_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/summarize_interactions/rank_across_celltypes_by_pair_by_avgsigmean/20200925.v3/cellphonedb.summary_across_celltypes_by_pair.min5.tsv")
## input druggable genes
drug_genes_df <- readxl::read_excel(path = "./Resources/Knowledge/Gene_Lists/Targetable_Genes.20200924.xlsx")

# filter by gene symbols --------------------------------------------------
cellphone_filtered_df <- cellphone_sum_by_pair_df %>%
  filter(rank_byinteraction_acrosscelltypes == 1) %>%
  mutate(gene.drug = ifelse(gene.source %in% drug_genes_df$genesymbol, gene.source,
                            ifelse(gene.target %in% drug_genes_df$genesymbol, gene.target, ""))) %>%
  filter(gene.drug != "")
  
cellphone_filtered_df$therapy_category <- mapvalues(x = cellphone_filtered_df$gene.drug, from = drug_genes_df$genesymbol, to = as.vector(drug_genes_df$therapy_category))
cellphone_filtered_df$pathway <- mapvalues(x = cellphone_filtered_df$gene.drug, from = drug_genes_df$genesymbol, to = as.vector(drug_genes_df$pathway))

# filter out duplicates ---------------------------------------------------



# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "filtered_druggale_related_pair_celltypes.", "tsv")
write.table(x = cellphone_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
