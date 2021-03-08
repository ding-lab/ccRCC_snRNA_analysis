# Yige Wu @WashU Mar 2021

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

# input denpendencies -----------------------------------------------------
## input marker table
markers1_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_germline/findallmarker_wilcox_germline_vhl_vs_1_somatic_nonsense_bycelltype_on_katmai/20200604.v1/VHL_Germline_vs_VHL_Somatic.FindMarkers.Wilcox.20200604.v1.tsv")
markers2_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_germline/findallmarker_wilcox_germline_vhl_vs_2nd_somatic_nonsense_bycelltype_on_katmai/20200604.v1/VHL_Germline_vs_VHL_Somatic.FindMarkers.Wilcox.20200604.v1.tsv")
markers2_df$group2_aliquot <- "CPT0015810004"
markers3_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_germline/findallmarker_wilcox_germline_vhl_vs_3rd_somatic_nonsense_bycelltype_on_katmai/20200604.v1/VHL_Germline_vs_VHL_Somatic.FindMarkers.Wilcox.20200604.v1.tsv")
markers3_df$group2_aliquot <- "CPT0086350004"
markers4_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_germline/findallmarker_wilcox_germline_vhl_vs_4th_somatic_nonsense_bycelltype_on_katmai/20200604.v1/VHL_Germline_vs_VHL_Somatic.FindMarkers.Wilcox.20200604.v1.tsv")
markers4_df$group2_aliquot <- "CPT0086820004"

# unite and transform markers -----------------------------------------------------------
markers_long_df <- rbind(markers1_df %>%
                           filter(Cell_type.shorter == "Tumor cells"),
                         markers2_df %>%
                           filter(Cell_type.shorter == "Tumor cells"))
markers_long_df <- rbind(markers_long_df,
                         markers3_df %>%
                           filter(Cell_type.shorter == "Tumor cells"))
markers_long_df <- rbind(markers_long_df,
                         markers4_df %>%
                           filter(Cell_type.shorter == "Tumor cells"))
## tranform
markers_filtered_long_df <- markers_long_df %>%
  filter(p_val_adj < 0.05)
markers_filtered_wide_df <- dcast(data = markers_filtered_long_df, formula = deg_gene_symbol ~ group2_aliquot, value.var = "avg_logFC")
markers_filtered_wide_df$number_up_VHL_Germline <- rowSums(markers_filtered_wide_df[, unique(markers_long_df$group2_aliquot)] > 0, na.rm = T)
markers_filtered_wide_df$number_down_VHL_Germline <- rowSums(markers_filtered_wide_df[, unique(markers_long_df$group2_aliquot)] < 0, na.rm = T)
markers_filtered_wide_df %>%
  filter(number_up_VHL_Germline == 4) %>%
  nrow()
markers_filtered_wide_df %>%
  filter(number_down_VHL_Germline == 4) %>%
  nrow()

# annotate those with top fold changes ------------------------------------
markers_wide_consistent_df <- markers_filtered_wide_df %>%
  filter(number_up_VHL_Germline == 4 | number_down_VHL_Germline == 4) %>%
  mutate(over3rdQu = ((CPT0001220012 >= quantile(x = markers_filtered_wide_df$CPT0001220012, 0.75, na.rm = T)) & (CPT0015810004 >= quantile(x = markers_filtered_wide_df$CPT0015810004, 0.75, na.rm = T)) & (CPT0086350004 >= quantile(x = markers_filtered_wide_df$CPT0086350004, 0.75, na.rm = T))& (CPT0086820004 >= quantile(x = markers_filtered_wide_df$CPT0086820004, 0.75, na.rm = T)))) %>%
  mutate(below1stQu = ((CPT0001220012 <= quantile(x = markers_filtered_wide_df$CPT0001220012, 0.25, na.rm = T)) & (CPT0015810004 <= quantile(x = markers_filtered_wide_df$CPT0015810004, 0.25, na.rm = T)) & (CPT0086350004 <= quantile(x = markers_filtered_wide_df$CPT0086350004, 0.25, na.rm = T))& (CPT0086820004 <= quantile(x = markers_filtered_wide_df$CPT0086820004, 0.25, na.rm = T))))

markers_wide_consistent_df %>%
  filter(over3rdQu) %>%
  nrow()
markers_wide_consistent_df %>%
  filter(below1stQu) %>%
  nrow()

# annotate with pathway ---------------------------------------------------
ora_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/findmarkers_germline/clusterprofiler/clusterprofiler_C3L-00917-T1_vs_4_somatic_nonsense_VHL_mutants_tumorcells/20210303.v1/ORA.Result.20210303.v1.tsv")
genes_ora <- sapply(ora_df$geneID, function(gene_string) {
  genes_tmp <- str_split(string = gene_string, pattern = "\\/")[[1]]
  return(genes_tmp)
})
genes_ora_uniq <- unique(unlist(genes_ora))
markers_wide_consistent_df$is_gene_in_ora_pathways <- (markers_wide_consistent_df$deg_gene_symbol %in% genes_ora_uniq)

View(markers_wide_consistent_df[markers_wide_consistent_df$is_gene_in_ora_pathways & (markers_wide_consistent_df$over3rdQu | markers_wide_consistent_df$below1stQu),])

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Filtered_DEGs_Wide.", run_id, ".tsv")
write.table(x = markers_filtered_wide_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "Consistent_DEGs_Wide.", run_id, ".tsv")
write.table(x = markers_wide_consistent_df, file = file2write, quote = F, sep = "\t", row.names = F)
