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

# input dependencies ------------------------------------------------------
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/findallmarker_wilcox_each_atac_tumor_vs_pt_on_katmai/20201130.v2/findallmarkers_wilcox_each_snatac_tumor_vs_pt.20201130.v2.tsv")
sample_anno_df <- data.frame(easyid = c("C3L-00088-N", "C3N-01200-N",
                                        "C3L-00416-T2", "C3L-01313-T1", "C3N-01200-T1",
                                        "C3L-00610-T1", "C3N-00733-T1",
                                        "C3L-00088-T1", "C3L-00088-T2", "C3L-00448-T1", "C3L-00917-T1"),
                             sample_group = c(rep("NAT", 2),
                                              rep("BAP1- tumor", 3),
                                              rep("PBRM1- tumor", 2),
                                              rep("non-mutant tumor", 4)))

# summarize genes by occurance ---------------------------------------------
unique(deg_df$easyid_tumor)
deg_filtered_df <- deg_df %>%
  filter(easyid_tumor %in% sample_anno_df$easyid) %>%
  filter(p_val_adj < 0.05) %>%
  filter(avg_logFC < 0)
deg_wide_df <- dcast(data = deg_filtered_df, formula = genesymbol_deg~easyid_tumor, value.var = "avg_logFC")
# countup_bap1_vec <- rowSums(x = (deg_wide_df[,sample_anno_df$easyid[sample_anno_df$sample_group == "BAP1- tumor"]] < 0), na.rm = T)
# countup_pbrm1_vec <- rowSums(x = (deg_wide_df[,sample_anno_df$easyid[sample_anno_df$sample_group == "PBRM1- tumor"]] < 0), na.rm = T)
# countup_nonmutant_vec <- rowSums(x = (deg_wide_df[,sample_anno_df$easyid[sample_anno_df$sample_group == "non-mutant tumor"]] < 0), na.rm = T)
countup_vec <- rowSums(x = (deg_wide_df[,-1] < 0), na.rm = T)
table(countup_vec)
# deg_wide_filtered_df <- deg_wide_df[countup_vec >= 9,]
deg_wide_filtered_df <- deg_wide_df[countup_vec >= 5,]

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "DEGs_all_snatac_tumors_down.", run_id, ".tsv")
write.table(x = deg_wide_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)

