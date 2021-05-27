# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
library(viridis)
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input meta data
id_metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210423.v1/meta_data.20210423.v1.tsv")
## input methylation data
methy_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Methylation/Methylation_Subtype/methylation_subtype_3_features.sig.50.txt")
# methy_df <- fread(data.table = F, input = "~/Box/CPTAC_ccRCC/Data_Freeze_1.1/CPTAC_ccRCC_Discovery/DNA_methylation/")
## input methylylation probe group
probe_group_df1 <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Methylation/Methylation_Subtype/methylation.sig.50.probes.by.metSubtype_3.1.txt")
probe_group_df2 <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Methylation/Methylation_Subtype/methylation.sig.50.probes.by.metSubtype_3.2.txt")
probe_group_df3 <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Methylation/Methylation_Subtype/methylation.sig.50.probes.by.metSubtype_3.3.txt")

# preprocess --------------------------------------------------------------
probe_group_df <- rbind(probe_group_df1 %>%
                          mutate(Probe_subtype = "m1"),
                        probe_group_df2 %>%
                          mutate(Probe_subtype = "m2"))
probe_group_df <- rbind(probe_group_df,
                        probe_group_df3 %>%
                          mutate(Probe_subtype = "m3"))

# format matrix data --------------------------------------------------
## get dim names
probe_names <- colnames(methy_df)
probe_names <- probe_names[grepl(pattern = "cg", x = probe_names)]
plot_data_t_mat <- as.matrix(methy_df[,probe_names])
plot_data_mat <- t(plot_data_t_mat)
colnames(plot_data_mat) <- methy_df$V1
## make row label
rownames_plot <- rownames(plot_data_mat)
colnames_plot <- colnames(plot_data_mat)

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## colors for gene set scores
colors_heatmapbody = colorRamp2(seq(0, 1, 0.2), 
                                viridis(n = 6))
# colors for group
colors_ismut <- c("black", "white")
names(colors_ismut) <- c("TRUE", "FALSE")
color_gridline = "grey90"


# make row split ----------------------------------------------------------
row_split_vec <- mapvalues(x = rownames_plot, from = probe_group_df$Gene, to = as.vector(probe_group_df$Probe_subtype))

# make column split -------------------------------------------------------
column_split_vec <- mapvalues(x = colnames_plot, from = methy_df$V1, to = as.vector(methy_df$subtype))
# mapvalues(x = "C3L-00908-T", from = methy_df$V1, to = as.vector(methy_df$subtype))
# mapvalues(x = "C3L-00416-T", from = methy_df$V1, to = as.vector(methy_df$subtype))
# mapvalues(x = "C3L-00416-T", from = methy_df$V1, to = as.vector(methy_df$BAP1_somatic_mutation))
colnames_plot[(paste0(colnames_plot, 1) %in% id_metadata_df$Aliquot.snRNA.WU[id_metadata_df$snATAC_available]) & (colnames_plot %in% methy_df$V1[methy_df$BAP1_somatic_mutation != ""])]
colnames_plot[(paste0(colnames_plot, 1) %in% id_metadata_df$Aliquot.snRNA.WU[id_metadata_df$snATAC_available]) & (colnames_plot %in% methy_df$V1[methy_df$BAP1_somatic_mutation != ""]) & (colnames_plot %in% methy_df$V1[methy_df$subtype != "1"])]

# make column annotation --------------------------------------------------
## merge data
colanno_df <- methy_df[, c("VHL_somatic_mutation", "PBRM1_somatic_mutation", "BAP1_somatic_mutation")]
rownames(colanno_df) <- methy_df$V1
colanno_df <- colanno_df[colnames_plot,]
## make colors
colors_scores_list <- list()
for (colname_tmp in colnames(colanno_df)) {
  colanno_df[, colname_tmp] <- as.character(colanno_df[, colname_tmp] != "")
  colors_scores_list[[colname_tmp]] <- colors_ismut
}
colanno_df$with_snRNA = (gsub(x = colnames_plot, pattern = "\\-T", replacement = "") %in% id_metadata_df$Case[id_metadata_df$snATAC_available])
colanno_df$with_snATAC = (gsub(x = colnames_plot, pattern = "\\-T", replacement = "") %in% id_metadata_df$Case[id_metadata_df$snATAC_available])
colors_scores_list[["with_snRNA"]] <- colors_ismut
colors_scores_list[["with_snATAC"]] <- colors_ismut
## make mark annotatin

colanno_obj = HeatmapAnnotation(df = colanno_df, col = colors_scores_list,
                                foo = anno_mark(at = which(colanno_df$with_snATAC == "TRUE"), 
                                                labels = colnames_plot[colanno_df$with_snATAC == "TRUE"], labels_gp = gpar(fontsize = 8)),
                                annotation_name_gp = gpar(fontsize = 14), annotation_name_side = "left", show_legend = F)

# plot  ------------------------------------------------------------
p <- ComplexHeatmap::Heatmap(matrix = plot_data_mat, 
                             col = colors_heatmapbody,
                             na_col = color_na, #border = "black",
                             ## row
                             show_row_names = F, row_split = row_split_vec,
                             show_row_dend = F, 
                             ## column
                             show_column_dend = F, cluster_columns = T, column_split = column_split_vec,
                             top_annotation = colanno_obj,
                             show_column_names = F,
                             show_heatmap_legend = F)
p


# make legend list --------------------------------------------------------
list_lgd = list(
  Legend(labels = names(colors_ismut), labels_gp = gpar(fontsize = 14),
         title = "Mutated", title_gp = gpar(fontsize = 14),
         legend_gp = gpar(fill = colors_ismut), border = NA),
  Legend(col_fun = colors_heatmapbody, 
         title = "Gene set score", 
         title_gp = gpar(fontsize = 14),
         labels_gp = gpar(fontsize = 14),
         legend_width = unit(4, "cm"),
         legend_height = unit(4, "cm"),
         direction = "horizontal"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Methylation", ".png")
png(file2write, width = 1200, height = 1100, res = 150)
draw(object = p, 
     annotation_legend_side = "top", annotation_legend_list = list_lgd)
dev.off()
file2write <- paste0(dir_out, "Methylation", ".pdf")
pdf(file2write, width = 9, height = 10, useDingbats = F)
draw(object = p,
     annotation_legend_side = "top", annotation_legend_list = list_lgd)
dev.off()

