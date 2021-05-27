# Yige Wu @WashU May 2021

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
## input the chromvar data
chromvar_sgbap1_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/enriched_motifs/dam_analysis/extract_chromvar_assay/extract_chromvar_786O_sgBAP1/20210525.v1/786O_sgBAP1.chromvar.tsv")
chromvar_control_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/enriched_motifs/dam_analysis/extract_chromvar_assay/extract_chromvar_786O_control/20210525.v1/786O_control.chromvar.tsv")
## input the motifs
jaspar <- fread(data.table = F, input = './Resources/Knowledge/snATAC/JASPAR2020_motifs.txt')

# test by motif -----------------------------------------------------------
wilcoxon_stat_df=NULL
for (i in 1:nrow(chromvar_sgbap1_df)) {
  scores_group1 <- chromvar_sgbap1_df[i, -1]; scores_group1 <- as.numeric(scores_group1)
  scores_group2 <- chromvar_control_df[i, -1]; scores_group2 <- as.numeric(scores_group2)
  w_test=wilcox.test(scores_group1, scores_group2)
  ## calculate the mean score of each group
  mean_score1=mean(scores_group1,na.rm=TRUE)
  mean_score2=mean(scores_group2,na.rm=TRUE)
  ## combine results and store
  stat=cbind(chromvar_sgbap1_df$V1[i],w_test$p.value,mean_score1,mean_score2)
  wilcoxon_stat_df=rbind(wilcoxon_stat_df,stat)
}
wilcoxon_stat_df <- as.data.frame(wilcoxon_stat_df)
colnames(wilcoxon_stat_df) <- c("motif_id", "p_val", "mean_score.sgBAP1", "mean_score.control")
wilcoxon_stat_df$motif_name <- mapvalues(x = wilcoxon_stat_df$motif_id, from = jaspar$motif, to = as.vector(jaspar$motif.name))
wilcoxon_stat_df$fdr <- p.adjust(p = wilcoxon_stat_df$p_val, method = "fdr")
wilcoxon_stat_df$mean_score.sgBAP1 <- as.numeric(wilcoxon_stat_df$mean_score.sgBAP1)
wilcoxon_stat_df$mean_score.control <- as.numeric(wilcoxon_stat_df$mean_score.control)
wilcoxon_stat_df$diff_score <- (wilcoxon_stat_df$mean_score.sgBAP1 - wilcoxon_stat_df$mean_score.control)
wilcoxon_stat_df <- wilcoxon_stat_df %>%
  arrange(p_val)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "786O_sgBAP1_vs_Control.Diff_Motifs.", run_id, ".tsv")
write.table(x = wilcoxon_stat_df, file = file2write, quote = F, sep = "\t", row.names = F)
