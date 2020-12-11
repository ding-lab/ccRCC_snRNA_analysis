# Yige Wu @WashU Nov 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 3
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies --------------------------------------------------------------
## input the dam annotation
# dam_anno_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/enriched_motifs/annotate_motifs/annotate_motifs_by_snatactumorgroup_shared/20201201.v1/dam_each_tumor_vs_pt.annotated.tsv")
dam_anno_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/enriched_motifs/annotate_motifs/annotate_motifs_by_snatactumorgroup_shared/20201202.v1/dam_each_tumor_vs_pt.annotated.tsv")

# filter ------------------------------------------------------------------
table(dam_anno_df$category_byshared)
## specify genes
dam_anno_df <- dam_anno_df %>%
  mutate(category_byshared_label = ifelse(category_byshared == "TRUE_TRUE_TRUE", "all-tumor-shared",
                                          ifelse(category_byshared == "FALSE_TRUE_TRUE", "PBRM-BAP1-shared",
                                                 ifelse(category_byshared == "FALSE_FALSE_TRUE", "PBRM1-specific",
                                                        ifelse(category_byshared == "FALSE_TRUE_FALSE", "BAP1-specific",
                                                               ifelse(category_byshared == "TRUE_FALSE_FALSE", "non-mutant-specific", "Others"))))))
  
dam_anno_df$mean_diff_acrosstumorgroups <- rowMeans(x = dam_anno_df[, c("mean_diff.nonmutant","mean_diff.bap1mutant","mean_diff.pbrm1mutant")], na.rm = T)
dam_anno_df <- dam_anno_df %>%
  mutate(abs_mean_diff = abs(mean_diff_acrosstumorgroups))

dam_filtered_df <- dam_anno_df %>%
  filter(category_byshared %in% c("FALSE_TRUE_FALSE", "FALSE_FALSE_TRUE", "FALSE_TRUE_TRUE", "TRUE_TRUE_TRUE", "TRUE_FALSE_FALSE"))
  arrange(category_byshared)
dam_filtered_df <- dam_filtered_df[!duplicated(dam_filtered_df$TF_Name),]
dam_top_up_df <- dam_filtered_df %>%
  filter(direction_shared == "up") %>%
  filter(!(TF_Name %in% c("FOS::JUNB", "NFE2", "BACH2", "NFE2L1"))) %>%
  group_by(category_byshared_label) %>%
  top_n(n = 5, wt = abs_mean_diff)
dam_top_up_df <- as.data.frame(dam_top_up_df)
dam_top_up_df <- rbind(dam_top_up_df,
                       dam_anno_df %>%
                         filter(TF_Name %in% c("ARNT::HIF1A")))
dam_top_down_df <- dam_filtered_df %>%
  filter(direction_shared == "down") %>%
  filter(category_byshared_label == "all-tumor-shared") %>%
  top_n(n = 10, wt = abs_mean_diff)

dam_top_df <- rbind(dam_top_up_df, dam_top_down_df)

# write outpu -------------------------------------------------------------
file2write <- paste0(dir_out, "Top_DAMs_Tumorcells_vs_PT_ByTumorGroup.", run_id, ".tsv")
write.table(x = dam_top_df, file = file2write, sep = "\t", quote = F, row.names = F)

