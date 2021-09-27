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
peaks_anno_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/emt/map_motifs_to_peaks/map_motifs_to_EMT_vs_selectedEpithelial_da_peaks/20210927.v3/EMT_vs_selectedEpithelial_diff_peaks_to_motifs.20210927.v3.tsv")
## input degs
degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_subclusters/findmarkers_selected_EMTclusters_vs_epithelialclusters_katmai/20210924.v3/Selected_2EMTclusters_vs_5Epithelialclusters.logfc.threshold0.min.pct0.1.min.diff.pct0.AssaySCT.tsv")
## input degs
dam_df <- fread(data.table = F, input = "./Resources/snATAC_Processed_Data/Enriched_Motifs/EMT/Score_difference.EpithelialSelectedClusters_vs_Mesenchymal.20210924.tsv")
