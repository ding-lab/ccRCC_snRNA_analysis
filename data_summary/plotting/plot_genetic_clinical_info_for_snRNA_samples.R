# Yige Wu @WashU March 2020
## running on local

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input id meta data
id_metadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv", data.table = F)
## input te bulk genomics/methylation events
bulk_sn_omicsprofile_df <- fread(input = "./Resources/Analysis_Results/data_summary/merge_bulk_sn_profiles/20200512.v1/bulk_sn_omics_profile.20200512.v1.tsv", data.table = F)
## input clinical info
specimen_clinical_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/extract_specimen_clinical_data/20200717.v1/snRNA_ccRCC_Specimen_Clinicl_Data.20200717.v1.tsv")

# merge data --------------------------------------------------------------
## filter samples with either snRNA data
id_metadata_filtered_df <- id_metadata_df %>%
  filter(snRNA_available) %>%
  filter(Is_discovery_set) %>%
  filter(Sample_Type == "Tumor")
## merge id meta data with bulk omics data
plot_data_df <- merge(id_metadata_filtered_df, bulk_sn_omicsprofile_df, by.x = c("Case", "Aliquot.snRNA", "Aliquot.snRNA.WU"), by.y = c("Case", "Aliquot.snRNA", "Aliquot_snRNA_WU"), all.x = T)
## merge clinical info
plot_data_df <- merge(plot_data_df,
                      specimen_clinical_df %>%
                        select(Sample, Histologic_Grade, Histologic_Type), 
                      by = c("Sample"), all.x = T)
plot_data_df <- plot_data_df %>%
  mutate(Aliquot_Suffix = str_split_fixed(string = Aliquot.snRNA.WU, pattern = "-", n = 3)[,3])
# plot_data_df$Aliquot_Suffix
# plot_data_df$Aliquot.snRNA.WU[order(match(plot_data_df$Aliquot_Suffix, c("T1", "T2", "T3", "N")))]
plot_data_df$snATAC_available <- (plot_data_df$Case %in% id_metadata_df$Case[id_metadata_df$snATAC_available])

# make data matrix for heatmap body ---------------------------------------
## reformat data frame to matrix
aliquot_ids <- plot_data_df$Aliquot.snRNA
plot_data_mat <- matrix(nc = length(aliquot_ids), nr = 0)
plot_data_mat %>% head()
### get case name
case_ids <- plot_data_df$Case
aliquot_wu_ids <- plot_data_df$Aliquot.snRNA.WU


# make colors -------------------------------------------------------------
## make colors for data availability
colors_bulk_data <- c(RColorBrewer::brewer.pal(n = 8, name = "Set1")[5], "white"); names(colors_bulk_data) <- c("TRUE", "FALSE")
colors_sn_data <- c(RColorBrewer::brewer.pal(n = 8, name = "Set1")[4], "white"); names(colors_sn_data) <- c("TRUE", "FALSE")
## make colors for histologic grade
colors_hist_grade <- c("NAT" = "white", "G1" = "#ffffcc", "G2" = "#addd8e", "G3" = "#31a354", "G4" = "#006837")
## make colors for histogical type
colors_hist_type <- c("Normal Adjacent Tissue" = "#66c2a5", "Clear cell renal cell carcinoma" = "#fc8d62", "non-Clear cell renal cell carcinoma" = "#8da0cb")
## top column annotation object
color_gridline = "grey90"
### make color for tumor purity
tumorpurity_color_fun <-  colorRamp2(c(0, 0.5, 1),c("white", "yellow", "red"))
## make colors for cnv state
colors_cnvstate <- c("#e41a1c", "#377eb8", "grey80")
names(colors_cnvstate)  <- c("gain", "loss", "neutral")

# make top column annotation --------------------------------------------------
top_col_anno_df <- plot_data_df 
# rownames(top_col_anno_df) <- as.vector(plot_data_df$Aliquot.snRNA)
# top_col_anno_df <- top_col_anno_df[aliquot_ids,]
# rownames(top_col_anno_df) <- aliquot_ids
### make neutral variant info for the normal sample
normal_aliquot_ids <- id_metadata_df$Aliquot.snRNA[id_metadata_df$Sample_Type == "Normal"]
#### for volumns other than methylation, make
top_col_anno_df[top_col_anno_df$Aliquot.snRNA %in% normal_aliquot_ids, grepl(x = colnames(top_col_anno_df), pattern = "Mut.")] <- "None"
top_col_anno_df[top_col_anno_df$Aliquot.snRNA %in% normal_aliquot_ids, paste0("CN.bulk.", c('3p', '5q', "14q"))] <- "neutral"
top_col_anno_df[top_col_anno_df$Aliquot.snRNA %in% normal_aliquot_ids, c("Translocation.t35", "Translocation.t32", "Translocation.t3_other")] <- "None"
top_col_anno_df[top_col_anno_df$Aliquot.snRNA %in% normal_aliquot_ids, c("Methyl.VHL")] <- NA
top_col_anno_df[top_col_anno_df$Aliquot.snRNA %in% normal_aliquot_ids, c("TumorPurity.snRNA", "TumorPurity.bulk")] <- 0
## make aliquot id suffix
aliquot_id_suffix <- str_split_fixed(string = aliquot_wu_ids, pattern = "-",n = 3)[,3]
aliquot_id_suffix
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
                                 quantile(top_col_anno_df$Methyl.VHL, 0.9, na.rm=T)),
                               c("#018571", "white", "#a6611a"))
### make text for translocations
Translocation.t3_other_text <- top_col_anno_df$Translocation.t3_other
Translocation.t3_other_text[is.na(Translocation.t3_other_text)] <- ""
Translocation.t3_other_text[Translocation.t3_other_text == "None"] <- ""
top_col_anno = HeatmapAnnotation(snRNA_Availability = anno_simple(x = as.character(top_col_anno_df$snRNA_available),
                                                         simple_anno_size = unit(3, "mm"),
                                                         gp = gpar(col = color_gridline), col = colors_sn_data),
                                 snATAC_Availability = anno_simple(x = as.character(top_col_anno_df$snATAC_available),
                                                          simple_anno_size = unit(3, "mm"), 
                                                          gp = gpar(col = color_gridline),  col = colors_sn_data),
                                 Histologic_Type = anno_simple(x = top_col_anno_df$Histologic_Type,
                                                               gp = gpar(col = color_gridline), 
                                                               simple_anno_size = unit(3, "mm"), 
                                                               col = colors_hist_type),
                                 Histologic_Grade = anno_simple(x = top_col_anno_df$Histologic_Grade,
                                                                gp = gpar(col = color_gridline),  
                                                                simple_anno_size = unit(3, "mm"), 
                                                                col = colors_hist_grade),
                                 gap1 = anno_empty(border = F, height = unit(0.5, "mm")),
                                 VHL_Methylation = anno_simple(x = top_col_anno_df$Methyl.VHL, 
                                                          gp = gpar(col = color_gridline), 
                                                          simple_anno_size = unit(3, "mm"),
                                                          col = methyl_color_fun),
                                 VHL = anno_simple(x = ifelse(is.na(top_col_anno_df$Mut.VHL), NA, 
                                                                  ifelse(!(top_col_anno_df$Mut.VHL == "None" | top_col_anno_df$Mut.VHL == "Silent"), 
                                                                         ifelse(top_col_anno_df$Is_discovery_set, "Mutated (WES)", "Mutated (Mapped the Mutation of T1 to snRNA Reads)"), "None")),
                                                       gp = gpar(col = color_gridline), 
                                                       simple_anno_size = unit(3, "mm"),
                                                       col = c("Mutated (WES)" = "#e7298a", "Mutated (Mapped the Mutation of T1 to snRNA Reads)" = "#c994c7", "None" = "white")),
                                 PBRM1 = anno_simple(x = ifelse(is.na(top_col_anno_df$Mut.PBRM1), NA, 
                                                                    ifelse(!(top_col_anno_df$Mut.PBRM1 == "None" | top_col_anno_df$Mut.PBRM1 == "Silent"), 
                                                                           ifelse(top_col_anno_df$Is_discovery_set, "Mutated (WES)", "Mutated (Mapped the Mutation of T1 to snRNA Reads)"), "None")),
                                                         gp = gpar(col = color_gridline), 
                                                         simple_anno_size = unit(3, "mm"),
                                                         col = c("Mutated (WES)" = "#e7298a", "Mutated (Mapped the Mutation of T1 to snRNA Reads)" = "#c994c7", "None" = "white")),
                                 BAP1 = anno_simple(x = ifelse(is.na(top_col_anno_df$Mut.BAP1), NA, 
                                                                   ifelse(!(top_col_anno_df$Mut.BAP1 == "None" | top_col_anno_df$Mut.BAP1 == "Silent"), 
                                                                          ifelse(top_col_anno_df$Is_discovery_set, "Mutated (WES)", "Mutated (Mapped the Mutation of T1 to snRNA Reads)"), "None")),
                                                        gp = gpar(col = color_gridline), 
                                                        simple_anno_size = unit(3, "mm"),
                                                        col = c("Mutated (WES)" = "#e7298a", "Mutated (Mapped the Mutation of T1 to snRNA Reads)" = "#c994c7", "None" = "white")),
                                 SETD2 = anno_simple(x = as.character(!(top_col_anno_df$Mut.SETD2 == "None" | top_col_anno_df$Mut.SETD2 == "Silent")),
                                                         gp = gpar(col = color_gridline), 
                                                         simple_anno_size = unit(3, "mm"),
                                                         col = c("TRUE" = "#e7298a", "FALSE" = "white")),
                                 KDM5C = anno_simple(x = as.character(!(top_col_anno_df$Mut.KDM5C == "None" | top_col_anno_df$Mut.KDM5C == "Silent")),
                                                         gp = gpar(col = color_gridline), 
                                                         simple_anno_size = unit(3, "mm"),
                                                         col = c("TRUE" = "#e7298a", "FALSE" = "white")),
                                 PTEN = anno_simple(x = as.character(!(top_col_anno_df$Mut.PTEN == "None" | top_col_anno_df$Mut.PTEN == "Silent")),
                                                        gp = gpar(col = color_gridline), 
                                                        simple_anno_size = unit(3, "mm"),
                                                        col = c("TRUE" = "#e7298a", "FALSE" = "white")),
                                 TSC1 = anno_simple(x = as.character(!(top_col_anno_df$Mut.TSC1 == "None" | top_col_anno_df$Mut.TSC1 == "Silent")),
                                                        gp = gpar(col = color_gridline), 
                                                        simple_anno_size = unit(3, "mm"),
                                                        col = c("TRUE" = "#e7298a", "FALSE" = "white")),
                                 gap2 = anno_empty(border = F, height = unit(0.5, "mm")),
                                 CN.3p = anno_simple(x = top_col_anno_df$CN.bulk.3p,
                                                          gp = gpar(col = color_gridline), 
                                                          simple_anno_size = unit(3, "mm"), 
                                                          col = colors_cnvstate),
                                 # CN.sn.3p_loss.fraction = anno_simple(top_col_anno_df$CN.sn.3p_loss.fraction,
                                 #                          simple_anno_size = unit(2, "mm")),
                                 CN.5q = anno_simple(x = top_col_anno_df$CN.bulk.5q,
                                                          gp = gpar(col = color_gridline), 
                                                          simple_anno_size = unit(3, "mm"), 
                                                          col = colors_cnvstate),
                                 # CN.sn.5q_gain.fraction = anno_barplot(x = top_col_anno_df$CN.sn.5q_gain.fraction, 
                                 #                                       height = unit(5, "mm")),
                                 CN.14q = anno_simple(x = top_col_anno_df$CN.bulk.14q,
                                                           gp = gpar(col = color_gridline), 
                                                           simple_anno_size = unit(3, "mm"), 
                                                           col = colors_cnvstate),
                                 annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 9))


# make column split -------------------------------------------------------
col_split_vec <- case_ids
## sort cases
case_sorted_df <- plot_data_df %>%
  mutate(Is_Mut.VHL = (Mut.VHL == "None")) %>%
  mutate(Is_Mut.PBRM1 = (Mut.PBRM1 == "None")) %>%
  mutate(Is_Mut.BAP1 = (Mut.BAP1 == "None")) %>%
  mutate(Is_Mut.SETD2 = (Mut.SETD2 == "None")) %>%
  mutate(Is_Mut.KDM5C = (Mut.KDM5C == "None")) %>%
  arrange(Histologic_Type, Is_Mut.VHL, Is_Mut.PBRM1, Is_Mut.BAP1, Is_Mut.SETD2, Is_Mut.KDM5C, snATAC_available)

## make factor
col_split_factor <- factor(x = case_ids, levels = case_sorted_df$Case)

# plot heatmap body with white-yellow-red ------------------------------------------------------
## make heatmap
color_na <- "grey50"
p <- Heatmap(matrix = plot_data_mat,
             na_col = color_na,
             ## column
             column_split = col_split_factor,
             column_order = order(match(plot_data_df$Aliquot_Suffix, c("T1", "T2", "T3", "N"))),
             column_title_gp = gpar(fontsize = 9), column_title_rot = 90,
             # bottom_annotation = bottom_col_anno, 
             top_annotation = top_col_anno,
             column_gap = unit(x = 0, units = "mm"),
             ## row
             show_heatmap_legend = F)
p
## make legend for top annotation
annotation_lgd = list(
  Legend(labels = c("Clear-cell RCC", "Non-clear-cell RCC"), labels_gp = gpar(fontsize = 9),
         title = "Histologic Type", title_gp = gpar(fontsize = 9),
         legend_gp = gpar(fill = colors_hist_type[-1]), border = color_gridline),
  Legend(labels = names(colors_hist_grade)[-1], labels_gp = gpar(fontsize = 9),
         title = "Histologic Grade (T1)", nrow = 2, title_gp = gpar(fontsize = 9),
         legend_gp = gpar(fill = colors_hist_grade[-1]), border = color_gridline),
  Legend(col_fun = methyl_color_fun, labels_gp = gpar(fontsize = 9),
         title = "VHL Promoter Methylation", title_gp = gpar(fontsize = 9),
         direction = "horizontal",
         legend_width = unit(30, "mm")),
  Legend(labels = c("Mutated (WES)", 
                    "None"), labels_gp = gpar(fontsize = 9),
         title = "Somatic Mutation Status", title_gp = gpar(fontsize = 9),
         legend_gp = gpar(fill = c("#e7298a", "white")), border = color_gridline),
  Legend(labels = c("Gain", "Loss", "Neutral"), labels_gp = gpar(fontsize = 9),
         title = "Bulk WGS CNV", nrow = 2, title_gp = gpar(fontsize = 9),
         legend_gp = gpar(fill = colors_cnvstate), border = color_gridline))
# ## save heatmap as png
# png(filename = paste0(dir_out, "data_availability",".png"), 
#     width = 1600, height = 1600, res = 150)
# ### combine heatmap and heatmap legend
# draw(object = p, 
#      annotation_legend_side = "right", annotation_legend_list = annotation_lgd)
# dev.off()
## save heatmap as pdf
pdf(paste0(dir_out, "data_availability", ".pdf"), 
    width = 6, height = 3.5)
### combine heatmap and heatmap legend
draw(object = p, 
     annotation_legend_side = "right", annotation_legend_list = annotation_lgd)
dev.off()

