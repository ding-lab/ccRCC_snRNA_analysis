# Yige Wu @WashU Jul 2020
## reference https://www.r-bloggers.com/2016/07/lets-be-faster-and-more-parallel-in-r-with-doparallel-package/

# set up libraries and output directory -----------------------------------
## getting the path to the current script
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
path_this_script <- thisFile()
## set working directory
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
## library additional libaries
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
library(doParallel)
## set run id
no_cores <- detectCores() - 1  
no_cores <- 25
# version_tmp <- "maxCores_minus1"
version_tmp <- paste0(no_cores, "Cores")
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
exp_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Methylation/CPTAC_ccRCC_discovery_tumor_methylation_betavalue_probe_level_v1.0.tsv")
## input bulk meta data
metadata_bulk_df <- fread("./Resources/Bulk_Processed_Data/Case_ID/CPTAC_ccRCC_discovery_caseID_v1.0.tsv")
## input mutation table
mut_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/mutation/annotate_cptac_sample_by_pbrm1_bap1_mutation/20210412.v1/PBRM1_BAP1_Mutation_Status_By_Case.20210412.v1.tsv")
## input methylation annotation
probe_anno_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Methylation/EPIC.hg38.manifest.tsv")

# get the aliquot IDs for bulk corresponding to the snRNA aliquots --------
metadata_filtered_df <- metadata_bulk_df %>%
  filter(Histologic_Type == "Clear cell renal cell carcinoma")
## add mutation group
metadata_filtered_df$group <- mapvalues(x = metadata_filtered_df$CASE_ID, from = mut_df$Case, to = as.vector(mut_df$mutation_category_sim))
metadata_filtered_df$group[as.vector(metadata_filtered_df$group) == as.vector(metadata_filtered_df$CASE_ID)] <- "Non-mutants"
table(metadata_filtered_df$group)
## get aliquout ids for the two groups
ids_exp_group1 <- metadata_filtered_df %>%
  filter(group %in% c("BAP1 mutated", "Both mutated")) %>%
  filter(!(CASE_ID %in% c("C3L-01287"))) %>%
  mutate(Sample.ID = paste0(CASE_ID, "-T"))
ids_exp_group1 <- ids_exp_group1$Sample.ID
ids_exp_group1 ## 16 samples
ids_exp_group2 <- metadata_filtered_df %>%
  filter(!(group %in% c("BAP1 mutated", "Both mutated"))) %>%
  mutate(Sample.ID = paste0(CASE_ID, "-T"))
ids_exp_group2 <- ids_exp_group2$Sample.ID
ids_exp_group2  ## 86 samples
## prepare probes to test
probe_test_df <- probe_anno_df %>%
  filter(!MASK_general) %>%
  filter(!is.na(gene_HGNC)) %>%
  filter(probeType == "cg") %>%
  filter(probeID %in% exp_df$Locus) %>%
  select(probeID)
nrow(probe_test_df)
probes_test <- head(probe_test_df$probeID, 5000)
probes_test <- probe_test_df$probeID
# test by wilcox and return values ----------------------------------------
exp_mat <- exp_df[,c(ids_exp_group1, ids_exp_group2)]
gene_index_vec = exp_df$Locus
ids_group1 = ids_exp_group1
ids_group2 = ids_exp_group2
# cl <- makeCluster(no_cores, type="FORK", outfile = paste0(dir_out, "makecluster.txt"))
# registerDoParallel(cl)
registerDoParallel(cores = no_cores)
start_time <- Sys.time()
test_list<-foreach(g=probes_test) %dopar% {
  exp_raw_vec1 <- unlist(exp_mat[gene_index_vec == g, ids_group1])
  exp_raw_vec2 <- unlist(exp_mat[gene_index_vec == g, ids_group2])
  exp_vec1 <- exp_raw_vec1[!is.na(exp_raw_vec1)]
  exp_vec2 <- exp_raw_vec2[!is.na(exp_raw_vec2)]
  if (length(exp_vec1) >= 5 & length(exp_vec2) >= 5) {
    log2FC <- log2(mean(exp_vec1)/mean(exp_vec2))
    stat <- wilcox.test(x = exp_vec1, y = exp_vec2)
    p_val <- stat$p.value
    result_list <- list(c(p_val, log2FC, length(exp_vec1), length(exp_vec2)))
  } else {
    result_list <- list(c(NA, NA, length(exp_vec1), length(exp_vec2)))
  }
  return(result_list)
}
end_time <- Sys.time()
end_time - start_time 
## 1.119908 mins for 1000 for 2 cores
## 37.37539 secs for 1000 for 4 cores
## 19.8362 for 1000 for 10 cores
## 51.0799 for 5000 for 20 cores
# stopCluster(cl)

## make test result into a data frame
test_df <- data.frame(matrix(data = unlist(test_list), ncol = 4, byrow = T))
colnames(test_df) <- c("p_val", "log2FC", "number_bap1tumors", "number_other_tumors")
test_df$probeID <- probes_test
## adjust p value
test_df$fdr <- p.adjust(p = test_df$p_val, method = "fdr")

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Methyaltion_BAP1_vs_Others.Wilcox.", run_id, ".tsv")
write.table(x = test_df, file = file2write, quote = F, sep = "\t", row.names = F)
