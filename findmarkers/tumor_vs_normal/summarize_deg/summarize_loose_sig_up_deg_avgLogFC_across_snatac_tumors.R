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
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/findallmarker_wilcox_each_atac_tumor_vs_pt_on_katmai/20201130.v2/findallmarkers_wilcox_each_snatac_tumor_vs_pt.20201130.v2.tsv")
sample_anno_df <- data.frame(easyid = c("C3L-00088-N", "C3N-01200-N",
                                        "C3L-01313-T1", "C3N-01200-T1", "C3L-01287-T1", "C3L-00416-T2", 
                                        "C3L-00610-T1", "C3N-00733-T1", "C3L-00079-T1", "C3L-00416-T2", 
                                        "C3L-00088-T1", "C3L-00088-T2", "C3L-00448-T1", "C3L-00917-T1"),
                             sample_group = c(rep("NAT", 2),
                                              rep("BAP1-mutant tumor", 4),
                                              rep("PBRM1-mutant tumor", 4),
                                              rep("non-mutant tumor", 4)))

# summarize genes by occurance ---------------------------------------------
unique(deg_df$easyid_tumor)
deg_filtered_df <- deg_df %>%
  filter(p_val_adj < 0.05) %>%
  filter(avg_logFC > 0)
deg_wide_df <- dcast(data = deg_filtered_df, formula = genesymbol_deg~easyid_tumor, value.var = "avg_logFC")

## summarize by tumor group
tumorgroup_tmp <- "BAP1-mutant tumor"
easyids_tmp <- as.vector(sample_anno_df$easyid[sample_anno_df$sample_group == tumorgroup_tmp])
deg_wide_df$count_bap1tumor_up <- rowSums(!is.na(deg_wide_df[, easyids_tmp]))
deg_wide_df$meanavglogFC_bap1tumor_up <- rowMeans(deg_wide_df[, easyids_tmp], na.rm = T)
tumorgroup_tmp <- "PBRM1-mutant tumor"
easyids_tmp <- as.vector(sample_anno_df$easyid[sample_anno_df$sample_group == tumorgroup_tmp])
deg_wide_df$count_pbrm1tumor_up <- rowSums(!is.na(deg_wide_df[, easyids_tmp]))
deg_wide_df$meanavglogFC_pbrm1tumor_up <- rowMeans(deg_wide_df[, easyids_tmp], na.rm = T)
tumorgroup_tmp <- "non-mutant tumor"
easyids_tmp <- as.vector(sample_anno_df$easyid[sample_anno_df$sample_group == tumorgroup_tmp])
deg_wide_df$count_nonmutanttumor_up <- rowSums(!is.na(deg_wide_df[, easyids_tmp]))
deg_wide_df$meanavglogFC_nonmutanttumor_up <- rowMeans(deg_wide_df[, easyids_tmp], na.rm = T)
## order
deg_wide_df <- deg_wide_df[, c("genesymbol_deg", 
                               "count_bap1tumor_up", "count_pbrm1tumor_up", "count_nonmutanttumor_up", 
                               "meanavglogFC_bap1tumor_up", "meanavglogFC_pbrm1tumor_up", "meanavglogFC_nonmutanttumor_up",
                               as.vector(sample_anno_df$easyid[sample_anno_df$sample_group != "NAT"]))]

# order -------------------------------------------------------------------
### order for bap1-specific
deg_wide_df <- deg_wide_df %>%
  arrange(desc(count_bap1tumor_up), count_pbrm1tumor_up, count_nonmutanttumor_up, desc(meanavglogFC_bap1tumor_up))
### order for all tumor common
deg_filtered_alltumorup_df <- deg_wide_df %>%
  filter(count_bap1tumor_up > 2 & count_pbrm1tumor_up > 2 & count_nonmutanttumor_up > 2) %>%
  arrange(desc(count_nonmutanttumor_up), desc(count_pbrm1tumor_up), desc(count_bap1tumor_up))
### order for BAP1 tumor specific
deg_filtered_bap1tumorup_df <- deg_wide_df %>%
  filter(count_bap1tumor_up > 2 & count_pbrm1tumor_up <=2 & count_nonmutanttumor_up <= 2) %>%
  arrange(desc(count_bap1tumor_up), count_pbrm1tumor_up, count_nonmutanttumor_up)
### order for PBRM1 tumor specific
deg_filtered_pbrm1tumorup_df <- deg_wide_df %>%
  filter(count_pbrm1tumor_up > 2 & count_bap1tumor_up <=2 & count_nonmutanttumor_up <=2) %>%
  arrange(desc(count_pbrm1tumor_up), count_bap1tumor_up, count_nonmutanttumor_up)
### order for BAP1/PBRM1 tumor shared
deg_filtered_bap1pbrm1tumorup_df <- deg_wide_df %>%
  filter(count_pbrm1tumor_up > 2 & count_bap1tumor_up > 2 & count_nonmutanttumor_up <=2) %>%
  arrange(desc(count_bap1tumor_up), desc(count_pbrm1tumor_up), count_nonmutanttumor_up)

# combine top motfs -------------------------------------------------------
deg_top_df <- rbind(deg_filtered_alltumorup_df %>%
                      mutate(DEG_category = "tumor-common-up"),
                    deg_filtered_bap1tumorup_df %>%
                      mutate(DEG_category = "BAP1-specific"))
deg_top_df <- rbind(deg_top_df,
                    deg_filtered_pbrm1tumorup_df %>%
                      mutate(DEG_category = "PBRM1-specific"))
deg_top_df <- rbind(deg_top_df,
                    deg_filtered_bap1pbrm1tumorup_df %>%
                      mutate(DEG_category = "BAP1/PBRM1-shared"))


# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "DEGs_sig_up.", run_id, ".tsv")
write.table(x = deg_wide_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "DEGs_top.", run_id, ".tsv")
write.table(x = deg_top_df, file = file2write, quote = F, sep = "\t", row.names = F)

