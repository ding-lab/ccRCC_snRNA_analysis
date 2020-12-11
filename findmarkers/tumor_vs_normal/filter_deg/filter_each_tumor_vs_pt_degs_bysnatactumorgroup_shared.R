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
deg_anno_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/annotate_deg/annotate_deg_by_snatactumorgroup_shared/20201202.v1/DEG_each_tumor_vs_pt.annotated.tsv")

# filter ------------------------------------------------------------------
deg_anno_df <- deg_anno_df %>%
  mutate(category_byshared_label = ifelse(category_byshared == "TRUE_TRUE_TRUE", "all-tumor-shared",
                                          ifelse(category_byshared == "FALSE_TRUE_TRUE", "PBRM-BAP1-shared",
                                                 ifelse(category_byshared == "FALSE_FALSE_TRUE", "PBRM1-specific",
                                                        ifelse(category_byshared == "FALSE_TRUE_FALSE", "BAP1-specific",
                                                               ifelse(category_byshared == "TRUE_FALSE_FALSE", "non-mutant-specific", "Others"))))))

deg_anno_df$mean_avg_logFC_acrosstumorgroups <- rowMeans(x = deg_anno_df[, c("mean_avg_logFC.nonmutant","mean_avg_logFC.bap1mutant","mean_avg_logFC.pbrm1mutant")], na.rm = T)
deg_anno_df <- deg_anno_df %>%
  mutate(abs_mean_avg_logFC = abs(mean_avg_logFC_acrosstumorgroups))
## filter up DEGs
deg_filtered_df <- deg_anno_df[!duplicated(deg_anno_df$genesymbol_deg),]
deg_top_up_df <- deg_filtered_df %>%
  filter(category_byshared_label %in% c("all-tumor-shared", "PBRM-BAP1-shared", "PBRM1-specific", "BAP1-specific", "non-mutant-specific")) %>%
  filter(direction_shared == "up") %>%
  filter(abs_mean_avg_logFC > 0.2) %>%
  # group_by(category_byshared_label) %>%
  # top_n(n = 50, wt = abs_mean_avg_logFC) %>%
  arrange(category_byshared_label, desc(abs_mean_avg_logFC))
deg_top_up_df <- as.data.frame(deg_top_up_df)
## filter down DEGs
deg_top_down_df <- deg_filtered_df %>%
  filter(category_byshared_label == "all-tumor-shared") %>%
  filter(direction_shared == "down") %>%
  filter(abs_mean_avg_logFC > 0.2) %>%
  # top_n(n = 50, wt = abs_mean_avg_logFC) %>%
  arrange(category_byshared_label, desc(abs_mean_avg_logFC))
deg_top_down_df <- as.data.frame(deg_top_down_df)
## combine
deg_top_df <- rbind(deg_top_up_df, deg_top_down_df)
nrow(deg_top_df)
table(deg_top_df$category_byshared_label)
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Top_DEGs_Tumorcells_vs_PT_ByTumorGroup.", run_id, ".tsv")
write.table(x = deg_top_df, file = file2write, quote = F, sep = "\t", row.names = F)

