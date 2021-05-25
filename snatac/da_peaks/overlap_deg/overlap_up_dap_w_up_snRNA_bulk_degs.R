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
deg_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/summarize_deg/overlap_tumor_vs_normal_snRNA_bulkRNA_protein_DEGs/20210511.v1/Tumor_vs_PT_DEGs.Overlap.snRNA.bulkRNA.Protein.20210511.v1.tsv")
dap_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/annotate_peaks/annotate_peaks_with_enhancer_prediction/20210511.v1/da_up_peaks_Tumor_vs_PT.annotated.wEnhancer.20210511.v1.tsv")

# overlap -----------------------------------------------------------------
deg_dap_merged_df <- merge(x = deg_df %>%
                             filter(direction.snRNA == "Up"), 
                           y = dap_df %>%
                             dplyr::arrange(desc(Count_up)) %>%
                             dplyr::rename(Count_up.snATAC = Count_up) %>%
                             dplyr::filter(Count_up.snATAC >= 5) %>%
                             dplyr::select(Gene, Count_up.snATAC, DAP_Type.strict, peak_distanceToTSS, peak), by.x = c("genesymbol_deg"), by.y = c("Gene"), all.x = T)
deg_dap_filtered_df <- deg_dap_merged_df %>%
  filter(DAP_Type.strict %in% c("Promoter", "Enhancer"))
deg_dap_promoter_df <- deg_dap_merged_df %>%
  filter(DAP_Type.strict == "Promoter")
deg_dap_enhancer_df <- deg_dap_merged_df %>%
  filter(DAP_Type.strict == "Enhancer")

# calculate ---------------------------------------------------------------
deg_df %>%
  filter(direction.snRNA == "Up") %>%
  nrow()
## 
table(deg_dap_filtered_df$genesymbol_deg)
unique(deg_dap_filtered_df$genesymbol_deg) %>% length()
##
deg_dap_promoter_df$genesymbol_deg %>% table()
unique(deg_dap_promoter_df$genesymbol_deg) %>% length()
deg_dap_enhancer_df$genesymbol_deg %>% table()
unique(deg_dap_enhancer_df$genesymbol_deg) %>% length()

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Up_DAP.Overlap.Up_snRNA_bulk_DEGs.", run_id, ".tsv")
write.table(x = deg_dap_merged_df, file = file2write, sep = "\t", row.names = F, quote = F)
file2write <- paste0(dir_out, "Up_Promoter_DAP.Overlap.Up_snRNA_bulk_DEGs.", run_id, ".tsv")
write.table(x = deg_dap_promoter_df, file = file2write, sep = "\t", row.names = F, quote = F)
file2write <- paste0(dir_out, "Up_ProEnh_DAP.Overlap.Up_snRNA_bulk_DEGs.", run_id, ".tsv")
write.table(x = deg_dap_filtered_df, file = file2write, sep = "\t", row.names = F, quote = F)
