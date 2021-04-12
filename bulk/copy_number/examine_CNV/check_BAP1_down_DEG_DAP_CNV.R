# Yige Wu @ WashU 2021 Ap
## annotate sample copy number profile (3p, 5q, 14q)

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/load_data.R")
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input CNA matrix
cna_df <- loadCNAstatus()
## input snRNA sample set
snRNA_meta_data <- fread("./Resources/Analysis_Results/sample_info/make_meta_data/20210322.v1/meta_data.20210322.v1.tsv", data.table = F)
## input the genes to inspect
dap_df <- fread(data.table = F, input = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/Differential_Peaks/BAP1_DOWN_top.20210408.witout1287.annotated.tsv")
## input mutation data
mut_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/mutation/annotate_cptac_sample_by_pbrm1_bap1_mutation/20210310.v1/PBRM1_BAP1_Mutation_Status_By_Case.20210310.v1.tsv")

# preprocess --------------------------------------------------------------
## specify the genes to be used 
genes2plot <- dap_df$Gene %>% unique()
## filter the CNVs
aliquots <- unique(snRNA_meta_data$Case[snRNA_meta_data$snRNA_available])
cna_filterd_wide_df <- cna_df[cna_df$gene %in% genes2plot, c("gene", aliquots)]
rownames(cna_filterd_wide_df) <- cna_filterd_wide_df$gene
cna_filterd_t_df <- t(cna_filterd_wide_df[, -1]) %>% as.data.frame()
cna_filterd_t_df$Case <- rownames(cna_filterd_t_df)
cna_filterd_t_df$mutation_category_sim <- mapvalues(x = cna_filterd_t_df$Case, from = mut_df$Case, to = as.vector(mut_df$mutation_category_sim))
cna_filterd_t_df <- cbind(cna_filterd_t_df[, c("Case", "mutation_category_sim")], cna_filterd_t_df[, colnames(cna_filterd_t_df)[!(colnames(cna_filterd_t_df) %in% c("Case", "mutation_category_sim"))]])
cna_filterd_t_df$mutation_category_sim[cna_filterd_t_df$Case == "C3L-00416"] <- "BAP1 mutated"
cna_filterd_t_df$mutation_category_sim[cna_filterd_t_df$Case == "C3L-01287"] <- "Non-mutants"

# count by gene by BAP1 status-----------------------------------------------------------
count_sup_df <- NULL
for (gene_tmp in cna_filterd_wide_df$gene) {
  count_df <- cna_filterd_t_df %>%
    mutate(BAP1_status = ifelse(mutation_category_sim == "BAP1 mutated", "BAP1_mutated", "BAP1_wt")) %>%
    select(BAP1_status, gene_tmp) %>%
    table() %>%
    as.data.frame() %>%
    mutate(genesymbol = gene_tmp)
  colnames(count_df)[colnames(count_df) == gene_tmp] <- "CNV_status"
  count_sup_df <- rbind(count_df, count_sup_df)
}
count_sup_df$peak <- mapvalues(x = count_sup_df$genesymbol, from = dap_df$Gene, to = as.vector(dap_df$peak))

# filter to the extreme ones ----------------------------------------------
cutoff_bap1mutant_num <- length(which(cna_filterd_t_df$mutation_category_sim == "BAP1 mutated"))/3
count_sup_filterd_df <- count_sup_df %>%
  filter(BAP1_status == "BAP1_mutated") %>%
  filter(CNV_status == "deletion") %>%
  filter(Freq > cutoff_bap1mutant_num) %>%
  arrange(desc(Freq))
count_sup_filterd_df$peak <- mapvalues(x = count_sup_filterd_df$genesymbol, from = dap_df$Gene, to = as.vector(dap_df$peak))
cutoff_bap1wt_num <- length(which(cna_filterd_t_df$mutation_category_sim != "BAP1 mutated"))/3
count_sup_filterd_df <- count_sup_df %>%
  filter(BAP1_status == "BAP1_wt") %>%
  filter(CNV_status == "amplification") %>%
  filter(Freq > cutoff_bap1wt_num) %>%
  arrange(desc(Freq))
count_sup_filterd_df$peak <- mapvalues(x = count_sup_filterd_df$genesymbol, from = dap_df$Gene, to = as.vector(dap_df$peak))

