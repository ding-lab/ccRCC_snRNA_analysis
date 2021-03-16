# Yige Wu @WashU Mar 2021
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
## input cellphonedb output
cellphone_df <- fread(data.table = F, input = "./Resources/Analysis_Results/cell_cell_interaction/annotate_interactions/annotate_cellphone_out/20201012.v1/cell.phone.res.total.run20200818.filtered.formatted.txt")
## input mutation category
mut_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/mutation/annotate_cptac_sample_by_pbrm1_bap1_mutation/20210310.v1/PBRM1_BAP1_Mutation_Status_By_Case.20210310.v1.tsv")

# map mutation status and summarize ---------------------------------------
## filter out normal sample
cellphone_filtered_df <- cellphone_df %>%
  filter(Sample_type == "Tumor") %>%
  mutate(Case = gsub(pattern = '\\-T[0-9]', replacement = "", x = Easy_id))
cellphone_filtered_df$mutation_category_sim <- mapvalues(x = cellphone_filtered_df$Case, from = mut_df$Case, to = as.vector(mut_df$mutation_category_sim))
sigcount_bypbrm1bap1_df <- cellphone_filtered_df %>%
  group_by(interacting_pair, pair_cell.types) %>%
  summarise(count_bap1mutant = sum(mutation_category_sim == "BAP1 mutated"),
            count_pbrm1mutant = sum(mutation_category_sim == "PBRM1 mutated"),
            count_nonmutant = sum(mutation_category_sim == "Non-mutants"),
            count_bothmutated = sum(mutation_category_sim == "Both mutated")) %>%
  as.data.frame()

# test --------------------------------------------------------------------
colnames_count <- c("count_bap1mutant", "count_pbrm1mutant", "count_nonmutant", "count_bothmutated")
names(colnames_count) <- c("BAP1 mutated", "PBRM1 mutated", "Non-mutants" ,"Both mutated")
samplecount <- cellphone_filtered_df %>%
  select(Easy_id, mutation_category_sim) %>%
  unique() %>%
  select(mutation_category_sim) %>%
  table()
group_tmp <- "BAP1 mutated"
for (group_tmp in c("BAP1 mutated", "PBRM1 mutated")) {
  pvalue_vec <- vector(mode = "numeric", length = nrow(sigcount_bypbrm1bap1_df))
  # i <- 1
  for (i in 1:nrow(sigcount_bypbrm1bap1_df)) {
    mat11 <- sigcount_bypbrm1bap1_df[i, colnames_count[group_tmp]]
    mat12 <- sum(sigcount_bypbrm1bap1_df[i, colnames_count[!(names(colnames_count) %in% group_tmp)]])
    mat21 <- samplecount[group_tmp] - mat11
    mat22 <- sum(samplecount) - samplecount[group_tmp] - mat12
    test_mat <- matrix(data = c(mat11, mat12, mat21, mat22), nrow = 2)
    fisher_out <- fisher.test(x = test_mat)
    pvalue_vec[i] <- fisher_out$p.value
  }
  sigcount_bypbrm1bap1_df[, gsub(pattern = "count", replacement = "pvalue", x = colnames_count[group_tmp])] <- pvalue_vec
  sigcount_bypbrm1bap1_df[, gsub(pattern = "count", replacement = "fdr", x = colnames_count[group_tmp])] <- p.adjust(p = pvalue_vec, method = "fdr")
}
sigcount_bypbrm1bap1_df <- sigcount_bypbrm1bap1_df %>%
  arrange(pvalue_bap1mutant, pvalue_pbrm1mutant)

# filter ------------------------------------------------------------------
sigcount_bypbrm1bap1_sig_df <- sigcount_bypbrm1bap1_df %>%
  filter(pvalue_bap1mutant < 0.05 | pvalue_pbrm1mutant < 0.05)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Fisher_Results.BAP1_PBRM1.", run_id, ".tsv")
write.table(x = sigcount_bypbrm1bap1_df,file = file2write, quote = F, sep = "\t", row.names = F)

