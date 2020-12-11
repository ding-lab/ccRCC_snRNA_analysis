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

# input dependencies --------------------------------------------------------------
## input dams
# dam_nonmutant_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/enriched_motifs/filter_motifs/filter_nonmutant_specific_motifs/20201201.v1/DAMs_nonmutant_snatac_tumors.tsv")
dam_nonmutant_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/enriched_motifs/filter_motifs/filter_nonmutant_specific_motifs/20201201.v2/DAMs_nonmutant_snatac_tumors.tsv")
dam_bap1_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/enriched_motifs/filter_motifs/filter_bap1_specific_motifs/20201201.v1/DAMs_bap1mutant_snatac_tumors.tsv")
dam_pbrm1_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/enriched_motifs/filter_motifs/filter_pbrm1_specific_motifs/20201201.v1/DAMs_pbrm1mutant_snatac_tumors.tsv")

# annotate ----------------------------------------------------------------
dam_anno_df <- dam_nonmutant_df %>%
  filter(!is.na(direction_shared)) %>%
  select(TF_Name, direction_shared)
dam_anno_df <- rbind(dam_anno_df,
                     dam_bap1_df %>%
                       filter(!is.na(direction_shared)) %>%
                       select(TF_Name, direction_shared))
dam_anno_df <- rbind(dam_anno_df,
                     dam_pbrm1_df %>%
                       filter(!is.na(direction_shared)) %>%
                       select(TF_Name, direction_shared))
dam_anno_df <- unique(dam_anno_df)
dam_anno_df <- merge(x = dam_anno_df, 
                     y = dam_nonmutant_df %>%
                       filter(!is.na(direction_shared)) %>%
                       mutate(mean_diff.nonmutant = mean_diff) %>%
                       select(TF_Name, direction_shared, mean_diff.nonmutant),
                     by = c("TF_Name", "direction_shared"), all.x = T)
dam_anno_df <- merge(x = dam_anno_df, 
                     y = dam_bap1_df %>%
                       filter(!is.na(direction_shared)) %>%
                       mutate(mean_diff.bap1mutant = mean_diff) %>%
                       select(TF_Name, direction_shared, mean_diff.bap1mutant),
                     by = c("TF_Name", "direction_shared"), all.x = T)
dam_anno_df <- merge(x = dam_anno_df, 
                     y = dam_pbrm1_df %>%
                       filter(!is.na(direction_shared)) %>%
                       mutate(mean_diff.pbrm1mutant = mean_diff) %>%
                       select(TF_Name, direction_shared, mean_diff.pbrm1mutant),
                     by = c("TF_Name", "direction_shared"), all.x = T)
dam_anno_df <- dam_anno_df %>%
  mutate(category_byshared = paste0(!is.na(mean_diff.nonmutant), "_", !is.na(mean_diff.bap1mutant), "_", !is.na(mean_diff.pbrm1mutant))) %>%
  arrange(desc(category_byshared), desc(mean_diff.nonmutant), desc(mean_diff.bap1mutant), desc(mean_diff.pbrm1mutant))

nrow(dam_anno_df)
length(unique(dam_anno_df$TF_Name))
table(dam_anno_df$category_byshared)
# FALSE_FALSE_TRUE FALSE_TRUE_FALSE  FALSE_TRUE_TRUE TRUE_FALSE_FALSE  TRUE_FALSE_TRUE  TRUE_TRUE_FALSE   TRUE_TRUE_TRUE 
# 641              766              719              142              213               61              481 

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "dam_each_tumor_vs_pt.annotated.", "tsv")
write.table(x = dam_anno_df, file = file2write, quote = F, sep = "\t", row.names = F)



