# Yige Wu @WashU Nov 2020

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
## input motif result
### mean score 1 is for PT, mean score 2 is for the tumor cells
dam_df <- fread(data.table = F, input = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/Enriched_Motifs/Tumor_vs_NormalPT/Score_difference.Tumor_Normal_comparison.20210322.tsv")
##
sample_anno_df <- data.frame(easyid = c("C3L-00088-N", "C3N-01200-N",
                                        "C3L-01313-T1", "C3N-01200-T1", "C3L-01287-T1", "C3L-00416-T2", "C3N-00317-T1",
                                        "C3L-00610-T1", "C3N-00733-T1", "C3L-00079-T1", "C3L-00416-T2", "C3N-00242-T1",
                                        "C3L-00088-T1", "C3L-00088-T2", "C3L-00448-T1", "C3L-00917-T1", "C3L-00096-T1"),
                             sample_group = c(rep("NAT", 2),
                                              rep("BAP1-mutant tumor", 5),
                                              rep("PBRM1-mutant tumor", 5),
                                              rep("non-mutant tumor", 5)))
## calculate cutoffs
sample_count_df <- sample_anno_df %>%
  select(sample_group) %>%
  table() %>%
  as.data.frame() %>%
  dplyr::rename(sample_group = '.')
cutoff_bap1 <- 0.5*sample_count_df$Freq[sample_count_df$sample_group == "BAP1-mutant tumor"]
cutoff_pbrm1 <- 0.5*sample_count_df$Freq[sample_count_df$sample_group == "PBRM1-mutant tumor"]
cutoff_nonmutant <- 0.5*sample_count_df$Freq[sample_count_df$sample_group == "non-mutant tumor"]

# summarize genes by occurance ---------------------------------------------
dam_filtered_df <- dam_df %>%
  mutate(easyid_tumor = gsub(x = cell_t2, pattern = "Tumor_", replacement = "")) %>%
  filter(FDR < 0.05) %>%
  filter(mean_score1 > 0 & diff < 0)
dam_wide_df <- dcast(data = dam_filtered_df, formula = TF_Name~easyid_tumor, value.var = "diff")
## summarize by tumor group
tumorgroup_tmp <- "BAP1-mutant tumor"
easyids_tmp <- as.vector(sample_anno_df$easyid[sample_anno_df$sample_group == tumorgroup_tmp])
dam_wide_df$count_bap1tumor_down <- rowSums(!is.na(dam_wide_df[, easyids_tmp]))
dam_wide_df$meandiff_bap1tumor_down <- rowMeans(dam_wide_df[, easyids_tmp], na.rm = T)
tumorgroup_tmp <- "PBRM1-mutant tumor"
easyids_tmp <- as.vector(sample_anno_df$easyid[sample_anno_df$sample_group == tumorgroup_tmp])
dam_wide_df$count_pbrm1tumor_down <- rowSums(!is.na(dam_wide_df[, easyids_tmp]))
dam_wide_df$meandiff_pbrm1tumor_down <- rowMeans(dam_wide_df[, easyids_tmp], na.rm = T)
tumorgroup_tmp <- "non-mutant tumor"
easyids_tmp <- as.vector(sample_anno_df$easyid[sample_anno_df$sample_group == tumorgroup_tmp])
dam_wide_df$count_nonmutanttumor_down <- rowSums(!is.na(dam_wide_df[, easyids_tmp]))
dam_wide_df$meandiff_nonmutanttumor_down <- rowMeans(dam_wide_df[, easyids_tmp], na.rm = T)
## order
dam_wide_df <- dam_wide_df[, c("TF_Name", 
                               "count_bap1tumor_down", "count_pbrm1tumor_down", "count_nonmutanttumor_down", 
                               "meandiff_bap1tumor_down", "meandiff_pbrm1tumor_down", "meandiff_nonmutanttumor_down",
                               as.vector(sample_anno_df$easyid[sample_anno_df$sample_group != "NAT"]))]


# order -------------------------------------------------------------------
### order for bap1-specific
dam_wide_df <- dam_wide_df %>%
  arrange(desc(count_bap1tumor_down), count_pbrm1tumor_down, count_nonmutanttumor_down)
### order for all tumor common
dam_filtered_alltumor_df <- dam_wide_df %>%
  filter(count_bap1tumor_down > cutoff_bap1 & count_pbrm1tumor_down > cutoff_pbrm1 & count_nonmutanttumor_down > cutoff_nonmutant)
dam_filtered_alltumor_df$order_score <- rowMeans(dam_filtered_alltumor_df[, sample_anno_df$easyid[sample_anno_df$sample_group != "NAT"]], na.rm = T)
dam_filtered_alltumor_df <- dam_filtered_alltumor_df %>%
  mutate(order_count = (count_nonmutanttumor_down+count_pbrm1tumor_down+count_bap1tumor_down)/3) %>%
  arrange(desc(order_count), desc(order_score))

# combine top motfs -------------------------------------------------------
dam_grouped_df <- dam_filtered_alltumor_df %>%
  mutate(TF_category = "tumor-common-down")

dam_top_df <- dam_grouped_df %>%
  head(n = 10)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "DAMs_sig_down.", run_id, ".tsv")
write.table(x = dam_wide_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "DAMs_sig_down.grouped.", run_id, ".tsv")
write.table(x = dam_grouped_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "DAMs_top_down.", run_id, ".tsv")
write.table(x = dam_top_df, file = file2write, quote = F, sep = "\t", row.names = F)
