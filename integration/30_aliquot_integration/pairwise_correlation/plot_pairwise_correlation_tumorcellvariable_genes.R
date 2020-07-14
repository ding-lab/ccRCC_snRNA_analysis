# Yige Wu @WashU May 2020
## running on local
## for plotting the aliquot-pairwise correlation coefficients for averaged expression of all genes

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
id_metadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200505.v1/meta_data.20200505.v1.tsv", data.table = F)
## input te bulk genomics/methylation events
bulk_sn_omicsprofile_df <- fread(input = "./Resources/Analysis_Results/bulk/other/merge_bulk_sn_profiles/20200512.v1/bulk_sn_omics_profile.20200512.v1.tsv", data.table = F)
## input the spearman pairwise correlation result
pearson_coef.tumorcellvariable_genes.df <- fread(input = "./Resources/Analysis_Results/integration/30_aliquot_integration/pairwise_correlation/calculate_pairwise_correlation_tumorcellvariable_genes/20200310.v1/avg_exp.tumorcellvaraible_genes.pearson_coef.20200310.v1.tsv", data.table = F)

# make data matrix for heatmap body ---------------------------------------
## reformat data frame to matrix
plot_data_df <- pearson_coef.tumorcellvariable_genes.df
plot_data_mat <- as.matrix(plot_data_df[,-1])
plot_data_mat %>% head()
rownames(plot_data_mat) <- plot_data_df$V1
### get case name
aliquot_ids <- rownames(plot_data_mat)
case_ids <- mapvalues(x = aliquot_ids, from = id_metadata_df$Aliquot.snRNA, to = as.vector(id_metadata_df$Case))
aliquot_wu_ids <- mapvalues(x = aliquot_ids, from = id_metadata_df$Aliquot.snRNA, to = as.vector(id_metadata_df$Aliquot.snRNA.WU))

# make top column annotation --------------------------------------------------
top_col_anno_df <- bulk_sn_omicsprofile_df %>%
  select(-Case) %>%
  select(-Aliquot.snRNA)
rownames(top_col_anno_df) <- as.vector(bulk_sn_omicsprofile_df$Aliquot.snRNA)
top_col_anno_df <- top_col_anno_df[aliquot_ids,]
rownames(top_col_anno_df) <- aliquot_ids
### make neutral variant info for the normal sample
normal_aliquot_ids <- id_metadata_df$Aliquot.snRNA[id_metadata_df$Sample_Type == "Normal"]
#### for volumns other than methylation, make
top_col_anno_df[rownames(top_col_anno_df) %in% normal_aliquot_ids, grepl(x = colnames(top_col_anno_df), pattern = "Mut.")] <- "None"
top_col_anno_df[rownames(top_col_anno_df) %in% normal_aliquot_ids, paste0("CN.", c('3p', '5q', "14q"))] <- "neutral"
top_col_anno_df[rownames(top_col_anno_df) %in% normal_aliquot_ids, c("Chr3_Translocation_Chr1", "Chr3_Translocation_Chr2")] <- "None"
top_col_anno_df[rownames(top_col_anno_df) %in% normal_aliquot_ids, c("Methyl.VHL")] <- NA
### make color for translocated chromosomes
uniq_translocation_chrs <- c(as.vector(top_col_anno_df$Chr3_Translocation_Chr1), as.vector(top_col_anno_df$Chr3_Translocation_Chr2))
uniq_translocation_chrs <- unique(uniq_translocation_chrs)
uniq_translocation_chrs
#### 7 different chromosomes
cartocolors_safe <- cartocolors_df[cartocolors_df$Name == "Safe", "n7"][[1]]
uniq_translocation_chr_colors <- c(cartocolors_safe, "white")
names(uniq_translocation_chr_colors) <- c(uniq_translocation_chrs[uniq_translocation_chrs!="None"], "None")
### make color for methylation
methyl_color_fun <- colorRamp2(c(quantile(top_col_anno_df$Methyl.VHL, 0.1, na.rm=T), 
                                 quantile(top_col_anno_df$Methyl.VHL, 0.5, na.rm=T), 
                                 quantile(top_col_anno_df$Methyl.VHL, 0.9, na.rm=T)),c("blue", "white", "red"))
### make color for tumor purity
tumorpurity_color_fun <-  colorRamp2(c(quantile(top_col_anno_df$TumorPurity.snRNA, 0.1, na.rm=T), 
                                       quantile(top_col_anno_df$TumorPurity.snRNA, 0.5, na.rm=T), 
                                       quantile(top_col_anno_df$TumorPurity.snRNA, 0.9, na.rm=T)),c("blue", "white", "red"))
### make text for translocations
Translocation.t3_other_text <- top_col_anno_df$Translocation.t3_other
Translocation.t3_other_text[is.na(Translocation.t3_other_text)] <- ""
Translocation.t3_other_text[Translocation.t3_other_text == "None"] <- ""
## top column annotation object
top_col_anno = HeatmapAnnotation(CN.bulk.3p = anno_simple(x = top_col_anno_df$CN.bulk.3p, 
                                                          gp = gpar(color = "black"),
                                                          simple_anno_size = unit(2, "mm"), 
                                                          col = cnv_state_colors),
                                 CN.sn.3p_loss.fraction = anno_barplot(top_col_anno_df$CN.sn.3p_loss.fraction,
                                                                       gp = gpar(fill = "lightblue", color = NA),
                                                                       height = unit(5, "mm")),
                                 CN.bulk.5q = anno_simple(x = top_col_anno_df$CN.bulk.5q,
                                                          gp = gpar(color = "black"),
                                                          simple_anno_size = unit(2, "mm"), 
                                                          col = cnv_state_colors),
                                 CN.sn.5q_gain.fraction = anno_barplot(x = top_col_anno_df$CN.sn.5q_gain.fraction, 
                                                                       gp = gpar(fill = "pink", color = NA),
                                                                       height = unit(5, "mm")),
                                 CN.bulk.14q = anno_simple(x = top_col_anno_df$CN.bulk.14q,
                                                           gp = gpar(color = "black"),
                                                           simple_anno_size = unit(2, "mm"), 
                                                           col = cnv_state_colors),
                                 CN.sn.14q_loss.fraction = anno_barplot(x = top_col_anno_df$CN.sn.14q_loss.fraction, 
                                                                        gp = gpar(fill = "lightblue", color = NA),
                                                                        height = unit(5, "mm")),
                                 Methyl.VHL = anno_simple(x = top_col_anno_df$Methyl.VHL, 
                                                          gp = gpar(color = "black"),
                                                          simple_anno_size = unit(3, "mm"),
                                                          col = methyl_color_fun),
                                 Mut.VHL = anno_simple(x = top_col_anno_df$Mut.VHL,
                                                       gp = gpar(color = "black"),
                                                       simple_anno_size = unit(3, "mm"),
                                                       col = variant_class_colors),
                                 Mut.PBRM1 = anno_simple(x = top_col_anno_df$Mut.PBRM1,
                                                         gp = gpar(color = "black"),
                                                         simple_anno_size = unit(3, "mm"),
                                                         col = variant_class_colors),
                                 Mut.SETD2 = anno_simple(x = top_col_anno_df$Mut.SETD2,
                                                         gp = gpar(color = "black"),
                                                         simple_anno_size = unit(3, "mm"),
                                                         col = variant_class_colors),
                                 Mut.BAP1 = anno_simple(x = top_col_anno_df$Mut.BAP1,
                                                        gp = gpar(color = "black"),
                                                        simple_anno_size = unit(3, "mm"),
                                                        col = variant_class_colors),
                                 Mut.KDM5C = anno_simple(x = top_col_anno_df$Mut.KDM5C,
                                                         gp = gpar(color = "black"),
                                                         simple_anno_size = unit(3, "mm"),
                                                         col = variant_class_colors),
                                 Mut.PTEN = anno_simple(x = top_col_anno_df$Mut.PTEN,
                                                        gp = gpar(color = "black"),
                                                        simple_anno_size = unit(3, "mm"),
                                                        col = variant_class_colors),
                                 Mut.TSC1 = anno_simple(x = top_col_anno_df$Mut.TSC1,
                                                        gp = gpar(color = "black"),
                                                        simple_anno_size = unit(3, "mm"),
                                                        col = variant_class_colors),
                                 Translocation.t35 = anno_simple(x = top_col_anno_df$Translocation.t35,
                                                                 gp = gpar(color = "black"),
                                                                 simple_anno_size = unit(3, "mm"),
                                                                 col = c("chr3-5" = "black", "None" = "white")),
                                 Translocation.t32 = anno_simple(x = top_col_anno_df$Translocation.t32,
                                                                 gp = gpar(color = "black"),
                                                                 simple_anno_size = unit(3, "mm"),
                                                                 col = c("chr3-2" = "black", "None" = "white")),
                                 Translocation.other = anno_text(x = Translocation.t3_other_text,
                                                                 which = "column", rot = 0, 
                                                                 location = 0.5, just = "center", height = unit(2, "mm"),
                                                                 gp = gpar(fill = "white", col = "black", fontsize = 5)),
                                 TumorPurity.snRNA = anno_simple(x = top_col_anno_df$TumorPurity.snRNA,
                                                                 gp = gpar(color = "black"),
                                                                 simple_anno_size = unit(3, "mm"),
                                                                 col = tumorpurity_color_fun))

# make color palette for each case ----------------------------------------
## input id meta data
uniq_case_ids <- unique(case_ids)
uniq_case_ids
### get unique color for each case
uniq_case_colors <- Polychrome::dark.colors(n = length(uniq_case_ids))
names(uniq_case_colors) <- uniq_case_ids

# make row annotation -------------------------------------------
## make case name as row annotation
### create row annotation
row_anno = rowAnnotation(foo = anno_text(aliquot_wu_ids, 
                                         location = 0.5, just = "center",
                                         # gp = gpar(fill = uniq_case_colors[case_ids], col = "white", border = "black"),
                                         gp = gpar(fill = NA, col = "black", border = NA),
                                         width = max_text_width(aliquot_wu_ids)*1.2))


# make bottom column annotation -------------------------------------------
bottom_col_anno = HeatmapAnnotation(foo = anno_text(aliquot_wu_ids, 
                                                    location = 0.5, just = "center",
                                                    # gp = gpar(fill = uniq_case_colors[case_ids], col = "white", border = "black"),
                                                    gp = gpar(fill = NA, col = "black", border = NA),
                                                    width = max_text_width(aliquot_wu_ids)*1.2))


# plot heatmap body with white-yellow-red ------------------------------------------------------
## make color function for heatmap body colors
col_fun = colorRamp2(c(0, 0.5, 1), c("white", "yellow", "red"))
## make heatmap
p <- Heatmap(matrix = plot_data_mat, 
             show_row_names = F, show_column_names = F, show_column_dend = F, show_row_dend = F,
             column_km = 2, column_km_repeats = 100,
             row_km = 2, row_km_repeats = 100,
             col = col_fun, 
             right_annotation = row_anno, 
             bottom_annotation = bottom_col_anno, 
             top_annotation = top_col_anno,
             show_heatmap_legend = F)
p
## make legend for heattmap body
heatmap_lgd = Legend(col_fun = col_fun, 
                     title = "Pearson's coeffcient\n(variably expressed genes\nwithin tumor cells)", 
                     direction = "vertical")
## make legend for top annotation
annotation_lgd = list(
  heatmap_lgd,
  Legend(labels = names(cnv_state_colors), 
         title = "Bulk WGS CNV", 
         legend_gp = gpar(fill = cnv_state_colors)),
  Legend(labels = names(variant_class_colors), 
         title = "Bulk Mutation Class", 
         legend_gp = gpar(fill = variant_class_colors)),
  Legend(labels = c("With Chr3 Translocation", "None"),
         title = "Bulk Chr3 Translocation",
         legend_gp = gpar(fill = c("black", "white"))),
  Legend(col_fun = methyl_color_fun,
         title = "Bulk VHL Promoter Methylation"),
  Legend(col_fun = tumorpurity_color_fun,
         title = "snRNA-based Tumor Purity"))
## save heatmap as png
png(filename = paste0(dir_out, "avg_exp.tumorcellvariable_genes.pearson_coef.heatmap.WYR.", run_id, ".png"), 
    width = 1200, height = 1800, res = 150)
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = annotation_lgd)
dev.off()

## save heatmap as pdf
pdf(file = paste0(dir_out, "avg_exp.tumorcellvariable_genes.pearson_coef.heatmap.WYR.", run_id, ".pdf"), 
    width = 8, height = 13)
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = annotation_lgd)
dev.off()


