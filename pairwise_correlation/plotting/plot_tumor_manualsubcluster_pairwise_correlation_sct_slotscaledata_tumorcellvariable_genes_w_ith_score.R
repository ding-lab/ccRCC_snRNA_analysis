# Yige Wu @WashU March 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input id meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv", data.table = F)
## input the spearman pairwise correlation result
pearson_coef.tumorcellvariable_genes.df <- fread(input = "./Resources/Analysis_Results/pairwise_correlation/calculate_tumor_manualsubcluster_pairwise_correlation_tumorcellvariable_genes_sct_scaledata/20201201.v1/avg_exp_by_tumorsubluster.tumorcellvaraible_genes.pearson_coef20201201.v1.tsv", data.table = F)
## input te bulk genomics/methylation events
bulk_sn_omicsprofile_df <- fread(input = "./Resources/Analysis_Results/data_summary/merge_bulk_sn_profiles/20200512.v1/bulk_sn_omics_profile.20200512.v1.tsv", data.table = F)
## input clinical info
specimen_clinical_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/extract_specimen_clinical_data/20200717.v1/snRNA_ccRCC_Specimen_Clinicl_Data.20200717.v1.tsv")
## input ITH score
ith_score_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_ith_score_assaysct_slotscaledata/20201201.v1/Mean_Pairwise_Corrleation_TumorSubclusters_By_Case.20201201.v1.tsv")

# make data matrix for heatmap body ---------------------------------------
## reformat data frame to matrix
plot_data_df <- pearson_coef.tumorcellvariable_genes.df
plot_data_mat <- as.matrix(plot_data_df[,-1])
dim(plot_data_mat)
plot_data_mat %>% head()
## add row names
rownames(plot_data_mat) <- plot_data_df$V1
plot_data_mat %>% head()

# process bulk omics ------------------------------------------------------
## filter samples with either snRNA data
id_metadata_filtered_df <- idmetadata_df %>%
  filter(snRNA_available)
## merge id meta data with bulk omics data
omics_df <- merge(id_metadata_filtered_df, bulk_sn_omicsprofile_df, by.x = c("Case", "Aliquot.snRNA", "Aliquot.snRNA.WU"), by.y = c("Case", "Aliquot.snRNA", "Aliquot_snRNA_WU"), all.x = T)
## merge clinical info
omics_df <- merge(omics_df,
                  specimen_clinical_df %>%
                    select(Sample, Histologic_Grade, Histologic_Type), 
                  by = c("Sample"), all.x = T)

rownames(omics_df) <- omics_df$Aliquot.snRNA.WU

# get ids -----------------------------------------------------------------
tumorsubcluster_ids <- rownames(plot_data_mat)
names_tumorsubclusters <- gsub(x = tumorsubcluster_ids, pattern = "\\.", replacement = "-")
names_tumorsubclusters
## get reable sample id
ids_aliquot_wu <- str_split_fixed(string = names_tumorsubclusters, pattern = "_", n = 2)[,1]
case_ids <- str_split_fixed(string = names_tumorsubclusters, pattern = "-T", n = 2)[,1]
## get suffixes of the aliquot ids
suffixes_aliquot_id <- str_split_fixed(string = ids_aliquot_wu, pattern = "-", n = 3)[,3]
suffixes_aliquot_id
names_cluster_suffix <- str_split_fixed(string = names_tumorsubclusters, pattern = "_", n = 2)[,2]
ids_showsuffix <- unique(case_ids[suffixes_aliquot_id == "T2"])
## sort case ids
### option 1: sort by number of subclusters
clustercount_by_case <- as.data.frame(table(case_ids))
clustercount_by_case <- clustercount_by_case %>%
  arrange(Freq)
# ids_case_bysubclustercount <- clustercount_by_case$case_ids
### option 2: sort by BAP1, PBRM1 and SETD2 mutation status
mutation_by_case <- omics_df %>%
  filter(Is_discovery_set) %>%
  filter(Sample_Type == "Tumor") %>%
  select(Case, Mut.VHL, Histologic_Type) %>%
  unique()
mutation_by_case$Mut.VHL[!is.na(mutation_by_case$Mut.VHL) & mutation_by_case$Mut.VHL != "None"] <- "Mutated"
# mutation_by_case$Mut.BAP1 <- factor(x = mutation_by_case$Mut.BAP1, levels = names(variant_class_colors))
# mutation_by_case$Mut.SETD2 <- factor(x = mutation_by_case$Mut.SETD2, levels = names(variant_class_colors))
# mutation_by_case$Mut.PBRM1 <- factor(x = mutation_by_case$Mut.PBRM1, levels = names(variant_class_colors))
# mutation_by_case <- mutation_by_case %>%
# arrange(Mut.BAP1, Mut.SETD2, Mut.PBRM1)
# ids_case_bymutation <- mutation_by_case$Case
### option 3: use both things
case_profile_df <- merge(clustercount_by_case, mutation_by_case, by.x = c("case_ids"), by.y = c("Case"))
case_profile_df$Mut.VHL <- factor(x = case_profile_df$Mut.VHL, levels = c("Mutated", "None"))
case_profile_df <- case_profile_df %>%
  arrange(Histologic_Type, Mut.VHL, Freq)
## order case ids
ids_case_uniq_ordered <- case_profile_df$case_ids
factor_case_ids <- factor(x = case_ids, levels = ids_case_uniq_ordered)
## get ITH scores
ith_scores_vec <- mapvalues(x = case_ids, from = ith_score_df$id_case1, to = as.vector(ith_score_df$ith_score))
ith_scores_vec <- as.numeric(ith_scores_vec)

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## make colors for histogical type
colors_hist_type <- c("Clear cell renal cell carcinoma" = "#fc8d62", "non-Clear cell renal cell carcinoma" = "#8da0cb")
## make color function for heatmap body colors
col_fun = colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red"))
## make color for the cluster name suffix
unique(names_cluster_suffix)
colors_clustername <- RColorBrewer::brewer.pal(n = 8, name = "Dark2")
names(colors_clustername) <- paste0("C", 1:8)
# swatch(colors_clustername)
## make colors for ITH score
summary(ith_scores_vec)
colors_ith_score <- colorRamp2(breaks = seq(0, 0.45, 0.05), colors = c("white", brewer.pal(n = 9, name = "YlGn")))

# make row annotation -------------------------------------------
## sort omics
omics_long_df <- omics_df[ids_aliquot_wu,]
## create left row annotation
row_anno_obj1 = rowAnnotation(foo = anno_block(gp = gpar(fill = "white", color = "black"),
                                               labels = ids_case_uniq_ordered, labels_rot = 0,
                                               labels_gp = gpar(col = "black", fontsize = 80)),
                              # Sample_Type_Suffix = anno_simple(x = suffixes_aliquot_id, col = colors_tumor_segments,
                              #                                  width = unit(2.5, "cm")),
                              annotation_name_side = "top", annotation_name_gp = gpar(fontsize = 20))


# make top column annotation -------------------------------------------
col_anno_obj1 = HeatmapAnnotation(Sample_Type_Suffix = anno_simple(x = suffixes_aliquot_id, col = colors_tumor_segments,
                                                                   height = unit(2.5, "cm")),
                                  ITH_Score = anno_simple(x = ith_scores_vec, col = colors_ith_score,
                                                          height = unit(2.5, "cm")),
                                  annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 20))

# plot pearson pairwise correlation for variably expressed genes within tumor cells ------------------------------------------------------
## save heatmap
p <- Heatmap(matrix = plot_data_mat,
             width = unit(nrow(plot_data_mat), "cm"), height = unit(ncol(plot_data_mat), "cm"),
             col = col_fun, 
             # use_raster = T, raster_device = "png",
             ## row
             row_split = factor_case_ids, cluster_row_slices = F,
             row_title = NULL,
             # row_title_rot = 0, row_title_gp = gpar(fontsize = 80),
             show_row_dend = F, 
             row_gap = unit(0, "mm"),
             # left_annotation = row_anno_obj2,
             left_annotation = row_anno_obj1,
             ## column
             column_split = factor_case_ids, cluster_column_slices = F, 
             column_title = NULL,
             # column_title_side = "bottom", column_title_rot = 90, column_title_gp = gpar(fontsize = 80),
             show_column_dend = F, 
             column_gap = unit(0, "mm"),
             # top_annotation = col_anno_obj2,
             top_annotation= col_anno_obj1,
             border = "grey50",
             show_row_names = F,
             show_column_names = F,
             show_heatmap_legend = F)
file2write <- paste0(dir_out, "Pariwise_Correlation", ".pdf")
pdf(file2write,
    width = 60, height = 60, useDingbats = F)
draw(p)
dev.off()


# plot legend -------------------------------------------------------------
## make horizontal legend
list_lgd = list(
  Legend(col_fun = col_fun, 
         title = "Pearson's coeffcient", 
         legend_width = unit(4, "cm"),
         direction = "horizontal"),
  Legend(col_fun = colors_ith_score, 
         title = "Intratumor heterogeneity\nscore by patient", 
         legend_width = unit(4, "cm"),
         direction = "horizontal"),
  Legend(title = "Tumor segment No.",
         labels = c("Tumor#1(T1)", "Tumor#2(T2)", "Tumor#3(T3)"),
         legend_gp = gpar(fill = colors_tumor_segments[c("T1", "T2", "T3")])),
  Legend(title = "Tumor Subcluster No.",
         labels = names(colors_clustername),
         legend_gp = gpar(fill = colors_clustername),
         nr = 2),
  Legend(labels = names(colors_hist_type),
         title = "Histologic Type",
         legend_gp = gpar(fill = colors_hist_type)))
dir_out_legend <- paste0(dir_out, "Legends/")
dir.create(dir_out_legend)
for (i in 1:length(list_lgd)) {
  file2write <- paste0(dir_out_legend, "Legend", i, ".pdf")
  pdf(file2write,
      width = 4, height = 3)
  draw(list_lgd[[i]])
  dev.off()
}

