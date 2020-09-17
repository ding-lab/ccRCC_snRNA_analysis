# Yige Wu @WashU Sep 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input tumor vs normal DA peaks
peaks2celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/da_peaks/filter_peaks/filter_celltype4_da_peaks_for_cellgoup4/20200916.v1/Cell_Type_Specific_DEGs_in_DA_peaks.20200916.v1.tsv")
## input peak-motif matrix
peaks2motif_df <- fread(data.table = F, input = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/4.Cell_type_markers/2.Update_annotation/Motif_matrix_MergedOb.tsv")
## input cell type specific motifs
motifs2celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/snatac/enriched_motifs/unite_motifs/unite_cellgroup4_enriched_motifs/20200915.v1/Enriched_Motifs.chromvar.MergedObj.byCell_group4.20200915.v1.tsv")

# filter cell type specifc peaks ------------------------------------------
peaks2celltype_filtered_df <- peaks2celltype_df %>%
  # mutate(pct_diff = (pct.1 - pct.2)) %>%
  # filter(pct_diff > 0.1) %>%
  group_by(Cell_type.filename) %>%
  mutate(rank_peak_by_avglogFC = order(order(avg_logFC.peak, decreasing = T)))

# filter cell type specific motifs ----------------------------------------
motifs2celltype_filtered_df <- motifs2celltype_df %>%
  filter(p_val_adj < 0.05) %>%
  mutate(pct_diff = (pct.1 - pct.2)) %>%
  group_by(Cell_type.filename) %>%
  mutate(rank_motif_by_avglogFC = order(order(abs(avg_logFC), decreasing = T))) %>%
  top_n(wt = abs(avg_logFC), n = 50) %>%
  select(-V1)

# motifs2celltype_filtered_df <- motifs2celltype_df %>%
#   filter(p_val_adj < 0.05) %>%
#   filter(pct_diff > 0.1) %>%
#   filter(Cell_type.filename %in% c("Tumor cells", "Normal epithelial cells")) %>%
#   group_by(Cell_type.filename) %>%
#   mutate(rank_motif_by_avglogFC = order(order(avg_logFC, decreasing = T))) %>%
#   top_n(wt = avg_logFC, n = 50) %>%
#   select(-V1)

# filter peak-motif matrix by pre-filtered peaks --------------------------
peaks2motif_df <- peaks2motif_df %>%
  filter(V1 %in% peaks2celltype_filtered_df$V1)
peaks2motif_melt_df <- melt(data = peaks2motif_df)
peaks2motif_filtered_df <- peaks2motif_melt_df %>%
  filter(value == 1) %>%
  filter(variable %in% motifs2celltype_filtered_df$motif) %>%
  select(-value)

# merge -------------------------------------------------------------------
peaks_merged_df <- merge(peaks2celltype_filtered_df, peaks2motif_filtered_df, by = c("V1"))
peaks_merged_df <- merge(peaks_merged_df, motifs2celltype_filtered_df, by.x = c("variable", "Cell_type.filename"), by.y = c("motif", "Cell_type.filename"), 
                         suffixes = c(".peak", ".motif"))
peaks_merged_df <- peaks_merged_df %>%
  arrange(Cell_type.filename, rank_peak_by_avglogFC, rank_motif_by_avglogFC)
unique(peaks_merged_df$SYMBOL)
unique(peaks2celltype_filtered_df$SYMBOL)
# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Cell_group4.DA_peaks..MotifAnnotated.", run_id, ".tsv")
write.table(x = peaks_merged_df, file = file2write, sep = "\t", quote = F, row.names = F)
