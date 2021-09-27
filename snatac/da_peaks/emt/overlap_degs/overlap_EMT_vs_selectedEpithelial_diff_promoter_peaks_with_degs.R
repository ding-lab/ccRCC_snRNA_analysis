# Yige Wu @WashU Sep 2021
## BAP1_tumorcells_vs_other_tumorcells

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
## input daps
peaks_anno_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/emt/annotate_peaks/annotate_EMT_vs_selectedEpithelial_diff_peaks_to_genes/20210927.v1/ccRCC_vs_PT_DAPs.Annotated.20210927.v1.tsv")
## input degs
degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/findmarkers_selected_EMTclusters_vs_epithelialclusters_katmai/20210924.v3/Selected_2EMTclusters_vs_5Epithelialclusters.logfc.threshold0.min.pct0.1.min.diff.pct0.AssaySCT.tsv")

# annotate and filter peaks ------------------------------------------------------------
## merge
peaks2degs_df <- merge(x = peaks_anno_df %>%
                         filter(Type == "Promoter") %>%
                         mutate(avg_log2FC.Epithelial_vs_Mesenchymal = avg_log2FC) %>%
                         mutate(avg_log2FC = -(avg_log2FC.Epithelial_vs_Mesenchymal)) %>%
                         select(peak, avg_log2FC, p_val_adj, Gene, Type, peak_distanceToTSS, pct.1, pct.2), 
                       y = degs_df,
                           by.x = c("Gene"), by.y = c("genesymbol_deg"), suffix = c(".snATAC", ".snRNA"), all = T)
peaks2degs_filtered_df <- peaks2degs_df %>%
  dplyr::filter(!is.na(avg_log2FC.snATAC) & !is.na(avg_log2FC.snRNA))

# write -------------------------------------------------------------------
file2write <- paste0(dir_out, "EMT_vs_selectedEpithelialClusters.", run_id, ".tsv")
write.table(x = peaks2degs_df, file = file2write, quote = F, sep = "\t", row.names = F)
