# Yige Wu @WashU Apr 2021

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
degs_ccRCC_vs_pt_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/summarize_deg/unite_tumor_vs_normal_snRNA_bulkRNA_protein_DEGs/20210608.v1/Tumor_vs_PT_DEGs.United.snRNA.bulkRNA.Protein.20210608.v1.tsv")
## input tumor-specific DEGs with surface annotation
# degs_ccRCC_vs_others_df
degs_ccRCC_vs_others_surface_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/tumor_specific_markers/20210701.v1/ccRCC_cells_specific_DEG_with_surface_annotations_from_3DB.txt")

# merge -------------------------------------------------------------------



# filter ------------------------------------------------------------------
degs_ccRCC_vs_otherspt_surface_df <- merge(x = degs_ccRCC_vs_others_surface_df, y = degs_ccRCC_vs_pt_df, by.x = "Gene", by.y = "genesymbol_deg", all.x = T)
degs_ccRCC_vs_otherspt_surface_filtered_df <- degs_ccRCC_vs_otherspt_surface_df %>%
  filter(GO_surface == "Surface") %>%
  filter(Num_sig_up >= 15 & Num_sig_down == 0)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "ccRCC_markers.Surface.", run_id, ".tsv")
write.table(x = degs_ccRCC_vs_otherspt_surface_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
