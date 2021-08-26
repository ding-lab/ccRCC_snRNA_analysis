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
degs_ccRCC_vs_pt_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/summarize_deg/unite_tumor_vs_normal_snRNA_bulkRNA_protein_DEGs/20210824.v1/Tumor_vs_PT_DEGs.United.snRNA.bulkRNA.Protein.20210824.v1.tsv")
## input tumor-specific DEGs with surface annotation
# degs_ccRCC_vs_others_df
degs_ccRCC_vs_others_surface_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/tumor_specific_markers/20210701.v1/ccRCC_cells_specific_DEG_with_surface_annotations_from_3DB.txt")
degs_ccRCC_vs_others_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/tumor_specific_markers/20210701.v1/Tumor cells_vs_combined_others_DE.txt")
## input druggable genes
druggenes_df1 <- fread(data.table = F, input = "./Resources/Knowledge/Gene_Lists/approved_target_ids_all.csv")
druggenes_df2 <- fread(data.table = F, input = "./Resources/Knowledge/Gene_Lists/01-Jul-2021-GeneSummaries.tsv")

# process selected genes -------------------------------------------------------------------
len_allsamples <- length(unique(degs_ccRCC_vs_others_df$sample_id))
degs_ccRCC_vs_others_sum_df <- degs_ccRCC_vs_others_df %>%
  group_by(gene_symbol) %>%
  summarise(avg_log2FC.mean.TumorcellsvsNontumor = mean(avg_log2FC), avg_norm_exp = mean(avg_norm_exp), sample_freq = length(sample_id)/len_allsamples) %>%
  rename(Gene = gene_symbol) %>%
  filter(Gene %in% c("CA9")) %>%
  mutate(GO_surface = NA) %>%
  mutate(CSPA_category = NA) %>%
  mutate(HPA_Reliability = NA)

# filter ------------------------------------------------------------------
markers_merged_df <- rbind(degs_ccRCC_vs_others_sum_df,
                                        degs_ccRCC_vs_others_surface_df %>%
                                          rename(avg_log2FC.mean.TumorcellsvsNontumor = avg_log2FC) %>%
                                          select(Gene, avg_log2FC.mean.TumorcellsvsNontumor, avg_norm_exp, sample_freq, GO_surface, CSPA_category, HPA_Reliability) %>%
                                          filter(GO_surface == "Surface" | (!is.na(HPA_Reliability) & HPA_Reliability %in% c("Approved", "Enhanced", "Supported"))))

markers_merged_df <- merge(x = markers_merged_df, 
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
markers_merged_filtered_df <- markers_merged_df %>%
  filter(Num_sig_up.allTumorcellsvsPT >= 15 & Num_sig_down.allTumorcellsvsPT == 0)
  
# annotate with druggability ----------------------------------------------
markers_merged_filtered_df <- merge(x = markers_merged_filtered_df, y = druggenes_df1 %>%
                                      select(`Gene Name`, `Drug IDs`) %>%
                                      rename(Drug_IDs = `Drug IDs`), by.x = c("Gene"), by.y = c("Gene Name"), all.x = T)
markers_merged_filtered_df <- merge(x = markers_merged_filtered_df, y = druggenes_df2 %>%
                                      select(name, description, gene_civic_url), by.x = c("Gene"), by.y = c("name"), all.x = T)

markers_merged_filtered_df <- markers_merged_filtered_df %>%
  mutate(is_druggable = (!is.na(Drug_IDs) | !is.na(gene_civic_url)))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "ccRCC_markers.Surface.", run_id, ".tsv")
write.table(x = markers_merged_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "ccRCC_markers.Surface.", run_id, ".csv")
write.table(x = markers_merged_filtered_df, file = file2write, quote = F, sep = ",", row.names = F)
file2write <- paste0(dir_out, "ccRCC_vs_NonTumorcells.merge.ccRCC_vs_PT.", run_id, ".csv")
write.table(x = markers_merged_df, file = file2write, quote = F, sep = ",", row.names = F)

# degs_ccRCC_vs_others_surface_df %>%
#   filter(GO_surface == "Surface" | (!is.na(HPA_Reliability) & HPA_Reliability %in% c("Approved", "Enhanced", "Supported"))) %>%
#   nrow()


