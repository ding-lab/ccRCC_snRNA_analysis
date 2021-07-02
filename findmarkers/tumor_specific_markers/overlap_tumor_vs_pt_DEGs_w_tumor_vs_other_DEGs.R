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
degs_ccRCC_vs_otherspt_surface_df <- merge(x = degs_ccRCC_vs_others_surface_df %>%
                                             rename(avg_log2FC.mean.TumorcellsvsNontumor = avg_log2FC) %>%
                                             select(Gene, avg_log2FC.mean.TumorcellsvsNontumor, avg_norm_exp, sample_freq, GO_surface, CSPA_category, HPA_Reliability), 
                                           y = degs_ccRCC_vs_pt_df %>%
                                             rename(pct.allTumorcells = pct.1.allTumorcellsvsPT) %>%
                                             rename(pct.allPTcells = pct.2.allTumorcellsvsPT) %>%
                                             rename(Num_sig_up.allTumorcellsvsPT = Num_sig_up) %>%
                                             rename(Num_sig_down.allTumorcellsvsPT = Num_sig_down) %>%
                                             rename(log2FC.bulkRNA = logFC.bulkRNA) %>%
                                             rename(log2FC.bulkpro = meddiff_exp.bulkpro) %>%
                                             select(genesymbol_deg, avg_log2FC.allTumorcellsvsPT, pct.allTumorcells, pct.allPTcells,
                                                    Num_sig_up.allTumorcellsvsPT, Num_sig_down.allTumorcellsvsPT, log2FC.bulkRNA, FDR.bulkRNA, log2FC.bulkpro, FDR.bulkpro), 
                                           by.x = "Gene", by.y = "genesymbol_deg", all.x = T)
degs_ccRCC_vs_otherspt_surface_filtered_df <- degs_ccRCC_vs_otherspt_surface_df %>%
  filter(GO_surface == "Surface" | (!is.na(HPA_Reliability) & HPA_Reliability %in% c("Approved", "Enhanced", "Supported"))) %>%
  filter(Num_sig_up.allTumorcellsvsPT >= 15 & Num_sig_down.allTumorcellsvsPT == 0)
  

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "ccRCC_markers.Surface.", run_id, ".tsv")
write.table(x = degs_ccRCC_vs_otherspt_surface_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)

# degs_ccRCC_vs_others_surface_df %>%
#   filter(GO_surface == "Surface" | (!is.na(HPA_Reliability) & HPA_Reliability %in% c("Approved", "Enhanced", "Supported"))) %>%
#   nrow()
