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
scores_wide_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/calculate_msigdb_geneset_scores/20211011.v1/MSigDB.Hallmark.tsv")

# preprocess --------------------------------------------------------------
## group gene sets into modules
module1_df <- data.frame(geneset_name = c("HALLMARK_MITOTIC_SPINDLE", "HALLMARK_E2F_TARGETS", "HALLMARK_G2M_CHECKPOINT", "HALLMARK_DNA_REPAIR", "HALLMARK_MYC_TARGETS_V1"),
                         module_name = "Cell_cycle")
module2_df <- data.frame(geneset_name = c("HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_COMPLEMENT", "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_KRAS_SIGNALING_UP"),
                         module_name = "Immune")
# module3_df <- data.frame(geneset_name = c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", "HALLMARK_HYPOXIA", "HALLMARK_TNFA_SIGNALING_VIA_NFKB"),
module3_df <- data.frame(geneset_name = c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"),
                         module_name = "EMT")
module4_df <- data.frame(geneset_name = c("HALLMARK_MTORC1_SIGNALING"),
                         module_name = "mTOR")
modules_df <- rbind(module1_df, module2_df, module3_df, module4_df)

# assign tumor cluster to group -----------------------------------------------------
for (module_tmp in unique(modules_df$module_name)) {
  genesetnames_tmp <- modules_df$geneset_name[modules_df$module_name == module_tmp]
  colnames_tmp <- paste0(gsub(x = genesetnames_tmp, pattern = "HALLMARK_", replacement = ""), "_Score")
  scores_tmp_df <- scores_wide_df[, colnames_tmp]; 
  scores_tmp_df <- data.frame(scores_tmp_df); colnames(scores_tmp_df) <- colnames_tmp
  rownames(scores_tmp_df) <- scores_wide_df$cluster_name
  
  if (module_tmp == "mTOR") {
    cutoffs_tmp <- quantile(x = scores_tmp_df[,colnames_tmp], probs = 0.9)
  } else {
    cutoffs_tmp <- sapply(scores_tmp_df, function(col_x) quantile(col_x, 0.75))
    if (module_tmp == "EMT") {
      cutoffs_tmp[grepl(x = names(cutoffs_tmp), pattern = "EPITHELIAL_MESENCHYMAL_TRANSITION", ignore.case = T)] <- quantile(x = scores_tmp_df[,"EPITHELIAL_MESENCHYMAL_TRANSITION_Score"], probs = 0.9)
    }
  }
  cutoffs_df <- matrix(data = rep(cutoffs_tmp, nrow(scores_tmp_df)), nrow = nrow(scores_tmp_df), byrow = T); cutoffs_df <- as.data.frame(cutoffs_df)
  scores_pass_df <- data.frame(scores_tmp_df > cutoffs_df)
  rownames(scores_pass_df) <- scores_wide_df$cluster_name
  count_pass_tmp <- rowSums(x = scores_pass_df)
  clusters_pass_tmp <- names(count_pass_tmp[count_pass_tmp == ncol(scores_pass_df)])
  scores_wide_df[, module_tmp] <- F; scores_wide_df[scores_wide_df$cluster_name %in% clusters_pass_tmp, module_tmp] <- T
}

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "MsigDB_Hallmark.Top15GeneSets.4Module.Enrichment.tsv")
write.table(x = scores_wide_df, file = file2write, quote = F, sep = "\t", row.names = F)
