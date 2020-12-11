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
## input degs
deg_nonmutant_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/filter_deg/filter_each_tumor_vs_normal_degs_to_nonmutant_shared/20201202.v1/DEGs_nonmutant_snatac_tumors.tsv")
deg_bap1_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/filter_deg/filter_each_tumor_vs_normal_degs_to_bap1mutant_shared/20201130.v1/DEGs_bap1mutant_snatac_tumors.tsv")
deg_pbrm1_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/filter_deg/filter_each_tumor_vs_normal_degs_to_pbrm1mutant_shared/20201202.v1/DEGs_pbrm1mutant_snatac_tumors.tsv")

# annotate ----------------------------------------------------------------
genesymbols_deg <- unique(c(deg_nonmutant_df$genesymbol_deg[!is.na(deg_nonmutant_df$direction_shared)], 
                            deg_bap1_df$genesymbol_deg[!is.na(deg_bap1_df$direction_shared)],
                            deg_pbrm1_df$genesymbol_deg[!is.na(deg_pbrm1_df$direction_shared)]))
deg_anno_df <- deg_nonmutant_df %>%
  filter(!is.na(direction_shared)) %>%
  select(genesymbol_deg, direction_shared)
deg_anno_df <- rbind(deg_anno_df,
                     deg_bap1_df %>%
                       filter(!is.na(direction_shared)) %>%
                       select(genesymbol_deg, direction_shared))
deg_anno_df <- rbind(deg_anno_df,
                     deg_pbrm1_df %>%
                       filter(!is.na(direction_shared)) %>%
                       select(genesymbol_deg, direction_shared))
deg_anno_df <- unique(deg_anno_df)
deg_anno_df <- merge(x = deg_anno_df, 
                     y = deg_nonmutant_df %>%
                       filter(!is.na(direction_shared)) %>%
                       mutate(mean_avg_logFC.nonmutant = mean_avg_logFC) %>%
                       select(genesymbol_deg, direction_shared, mean_avg_logFC.nonmutant),
                     by = c("genesymbol_deg", "direction_shared"), all.x = T)
deg_anno_df <- merge(x = deg_anno_df, 
                     y = deg_bap1_df %>%
                       filter(!is.na(direction_shared)) %>%
                       mutate(mean_avg_logFC.bap1mutant = mean_avg_logFC) %>%
                       select(genesymbol_deg, direction_shared, mean_avg_logFC.bap1mutant),
                     by = c("genesymbol_deg", "direction_shared"), all.x = T)
deg_anno_df <- merge(x = deg_anno_df, 
                     y = deg_pbrm1_df %>%
                       filter(!is.na(direction_shared)) %>%
                       mutate(mean_avg_logFC.pbrm1mutant = mean_avg_logFC) %>%
                       select(genesymbol_deg, direction_shared, mean_avg_logFC.pbrm1mutant),
                     by = c("genesymbol_deg", "direction_shared"), all.x = T)
deg_anno_df <- deg_anno_df %>%
  mutate(category_byshared = paste0(!is.na(mean_avg_logFC.nonmutant), "_", !is.na(mean_avg_logFC.bap1mutant), "_", !is.na(mean_avg_logFC.pbrm1mutant))) %>%
  arrange(desc(category_byshared), desc(mean_avg_logFC.nonmutant), desc(mean_avg_logFC.bap1mutant), desc(mean_avg_logFC.pbrm1mutant))

nrow(deg_anno_df)
length(unique(deg_anno_df$genesymbol_deg))
table(deg_anno_df$category_byshared)
# FALSE_FALSE_TRUE FALSE_TRUE_FALSE  FALSE_TRUE_TRUE TRUE_FALSE_FALSE  TRUE_FALSE_TRUE  TRUE_TRUE_FALSE   TRUE_TRUE_TRUE 
# 641              766              719              142              213               61              481 

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "DEG_each_tumor_vs_pt.annotated.", "tsv")
write.table(x = deg_anno_df, file = file2write, quote = F, sep = "\t", row.names = F)



