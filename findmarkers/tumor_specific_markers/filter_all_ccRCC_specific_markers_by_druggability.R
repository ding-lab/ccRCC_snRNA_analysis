# Yige Wu @WashU July 2021

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
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
markers_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_specific_markers/overlap_tumor_specifc_markers_w_tumor_vs_NAT_DEGs/20220711.v1/ccRCC_cells_specific_DEG_with_surface_annotations_from_3DB.filtred.txt")
## input druggable genes
druggenes_df1 <- fread(data.table = F, input = "./Resources/Knowledge/Gene_Lists/approved_target_ids_all.csv")
druggenes_df2 <- fread(data.table = F, input = "./Resources/Knowledge/Gene_Lists/01-Jul-2021-GeneSummaries.tsv")
## merge DEGs
tumorvspt_degs_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/findmarker_LR_all_ccRCC_vs_pt_on_katmai/20210824.v1/LR.logfc.threshold0.min.pct0.1.min.diff.pct0.AssayRNA.tsv")

# filter ------------------------------------------------------------------
markers_merged_df <- merge(x = markers_df, y = druggenes_df1, by.x = c("Gene"), by.y = c("Gene Name"), all.x = T)
markers_merged_df <- merge(x = markers_merged_df, y = druggenes_df2, by.x = c("Gene"), by.y = c("name"), all.x = T)
markers_merged_df <- merge(x = markers_merged_df, y = tumorvspt_degs_df %>%
                             select(genesymbol_deg, avg_log2FC, p_val_adj) %>%
                             rename(avg_log2FC.Tumorcell_vs_PT = avg_log2FC) %>%
                             rename(FDR.Tumorcell_vs_PT = p_val_adj), by.x = c("Gene"), by.y = c("genesymbol_deg"), all.x = T)

markers_filtered_df <- markers_merged_df %>%
  filter(!is.na(`Drug IDs`) | !is.na(gene_civic_url)) %>%
  arrange(desc(avg_log2FC)) %>%
  unique()

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "ccRCC_markers.Druggable.Annotated.", run_id, ".tsv")
write.table(x = markers_merged_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "ccRCC_markers.Druggable.", run_id, ".tsv")
write.table(x = markers_filtered_df, file = file2write, quote = F, sep = "\t", row.names = F)
