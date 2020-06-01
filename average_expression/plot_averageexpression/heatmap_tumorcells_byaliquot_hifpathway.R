# Yige Wu @WashU Apr 2020
## plot heatmap the the average expression (not scaled) of HIF pathway members to compare tumor and normal epithelial cells

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

# input denpendencies -----------------------------------------------------
## input average expression by cell type by aliquot
avgexp_bycelltype_byaliquot_df <- fread(input = "./Resources/Analysis_Results/average_expression/averageexpression_bycelltypedetailed_byaliquot_on_katmai/20200413.v1/averageexpression_bycelltypedetailed.30_aliquot_integration.20200413.v1.tsv", data.table = F)
## input the hif pathway members
hiftargets_df <- fread(input = "./Resources/Analysis_Results/dependencies/write_hif_targets/20200302.v1/HIF_Target_Genes.20200302.v1.tsv", data.table = F)
## input barcode 2 cell type table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/map_celltype_to_barcode/map_celltype_to_all_cells/20200410.v1/30_aliquot_integration.barcode2celltype.20200410.v1.tsv", data.table = F)
## input id meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200413.v1/meta_data.20200413.v1.tsv", data.table = F)
## input te bulk genomics/methylation events
bulk_sn_omicsprofile_df <- fread(input = "./Resources/Analysis_Results/bulk/other/merge_bulk_sn_profiles/20200512.v1/bulk_sn_omics_profile.20200512.v1.tsv", data.table = F)

# set plot gene list ------------------------------------------------------
genes_plot <- c(hiftargets_df$target_genesymbol, "HIF1A", "EPAS1", "VHL")

# count number of cells per celltype per aliquot --------------------------
count_bycelltype_byaliquot_df <- data.frame(table(barcode2celltype_df[, c("Cell_type.detailed", "orig.ident")]))
unique(count_bycelltype_byaliquot_df$Cell_type.detailed)
count_bycelltype_byaliquot_df$celltype <- gsub(x = count_bycelltype_byaliquot_df$Cell_type.detailed, pattern = '[/ +]', replacement = ".")
count_bycelltype_byaliquot_df$celltype <- gsub(x = count_bycelltype_byaliquot_df$celltype, pattern = '\\-', replacement = ".")
count_bycelltype_byaliquot_df <- count_bycelltype_byaliquot_df %>%
  mutate(id_aliquot_celltype = paste0(orig.ident, "_", celltype))
count_bycelltype_byaliquot_filtered <- count_bycelltype_byaliquot_df %>%
  filter(celltype %in% c("Tumor.cells")) %>%
  filter(Freq >= 10)

# make matrix for heatmap body --------------------------------------------
plot_data_df <- avgexp_bycelltype_byaliquot_df
## format the column names to only aliquot id + cell type
data_col_names <- colnames(plot_data_df)[-1]
data_col_names
data_col_names.changed <- str_split_fixed(string = data_col_names, pattern = "\\.", n = 2)[,2]
data_col_names.changed
## rename the data frame
colnames(plot_data_df) <- c("gene", data_col_names.changed)
## get the column names for tumor cells and normal epithelial cells with enough cells
data_col_names.keep <- count_bycelltype_byaliquot_filtered$id_aliquot_celltype
# data_col_names.keep <- data_col_names[grepl(pattern = "Normal.epithelial.cells", x = data_col_names) | grepl(pattern = "Tumor.cells", x = data_col_names)]
data_col_names.keep
## get the gene names for HIF targets
data_row_names <- plot_data_df$gene
data_row_names.keep <- intersect(data_row_names, genes_plot)
data_row_names.keep
## filter cell group down to epithelial cells
## reformat data frame to matrix
plot_data_mat <- as.matrix(plot_data_df[,data_col_names.keep])
plot_data_mat %>% head()
## add row names
rownames(plot_data_mat) <- data_row_names
plot_data_mat %>% head()
## filter down to HIF target genes
plot_data_mat <- plot_data_mat[data_row_names.keep,]
### get aliquot ids and case ids
ids_aliquot_celltype <- data_col_names.keep
ids_aliquot <- str_split_fixed(string = ids_aliquot_celltype, pattern = "_", n = 2)[,1]
ids_aliquot_wu <- mapvalues(x = ids_aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
ids_celltype <- str_split_fixed(string = ids_aliquot_celltype, pattern = "_", n = 2)[,2]
# case_ids <- mapvalues(x = ids_aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Case))
## change coluumn ids
colnames(plot_data_mat) <- ids_aliquot

# make column split -------------------------------------------------------
idaliquot_vhl_germline <- bulk_sn_omicsprofile_df$Aliquot.snRNA[!is.na(bulk_sn_omicsprofile_df$Mut.VHL.Germline) & (bulk_sn_omicsprofile_df$Mut.VHL.Germline != "None")]
idaliquot_vhl_germline
idaliquot_3pdel_only <- bulk_sn_omicsprofile_df$Aliquot.snRNA[!is.na(bulk_sn_omicsprofile_df$Mut.VHL) & (bulk_sn_omicsprofile_df$Mut.VHL == "None")]
columnsplit_vec <- ifelse(ids_aliquot %in% idaliquot_vhl_germline, "VHL_Germline",
                          ifelse(ids_aliquot %in% idaliquot_3pdel_only, "3p_Loss", "VHL_Somatic\n+3p_Loss"))
columnsplit_factor <- factor(x = columnsplit_vec, levels = c("VHL_Germline", "VHL_Somatic\n+3p_Loss", "3p_Loss"))

# make top column annotation -----------------------------------------------------
top_col_anno_df <- bulk_sn_omicsprofile_df %>%
  select(-Case) %>%
  select(-Aliquot.snRNA)
rownames(top_col_anno_df) <- as.vector(bulk_sn_omicsprofile_df$Aliquot.snRNA)
top_col_anno_df <- top_col_anno_df[ids_aliquot,]
rownames(top_col_anno_df) <- ids_aliquot
### make neutral variant info for the normal sample
normal_ids_aliquot <- idmetadata_df$Aliquot.snRNA[idmetadata_df$Sample_Type == "Normal"]
#### for volumns other than methylation, make
top_col_anno_df[rownames(top_col_anno_df) %in% normal_ids_aliquot, grepl(x = colnames(top_col_anno_df), pattern = "Mut.")] <- "None"
top_col_anno_df[rownames(top_col_anno_df) %in% normal_ids_aliquot, paste0("CN.", c('3p', '5q', "14q"))] <- "neutral"
top_col_anno_df[rownames(top_col_anno_df) %in% normal_ids_aliquot, c("Chr3_Translocation_Chr1", "Chr3_Translocation_Chr2")] <- "None"
top_col_anno_df[rownames(top_col_anno_df) %in% normal_ids_aliquot, c("Methyl.VHL")] <- NA
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
                                                          simple_anno_size = unit(2, "mm"), 
                                                          col = cnv_state_colors),
                                 CN.sn.3p_loss.fraction = anno_barplot(top_col_anno_df$CN.sn.3p_loss.fraction, 
                                                                       height = unit(5, "mm")),
                                 Methyl.VHL = anno_simple(x = top_col_anno_df$Methyl.VHL,
                                                          simple_anno_size = unit(3, "mm"),
                                                          col = methyl_color_fun),
                                 Mut.VHL.Somatic = anno_simple(x = top_col_anno_df$Mut.VHL,
                                                       simple_anno_size = unit(4, "mm"),
                                                       col = variant_class_colors),
                                 Mut.VHL.Germline = anno_simple(x = top_col_anno_df$Mut.VHL.Germline,
                                                       simple_anno_size = unit(4, "mm"),
                                                       col = variant_class_colors),
                                 CN.bulk.14q = anno_simple(x = top_col_anno_df$CN.bulk.14q,
                                                           simple_anno_size = unit(2, "mm"), 
                                                           col = cnv_state_colors),
                                 CN.sn.14q_loss.fraction = anno_barplot(x = top_col_anno_df$CN.sn.14q_loss.fraction, 
                                                                        height = unit(5, "mm")),
                                 Translocation.t35 = anno_simple(x = top_col_anno_df$Translocation.t35,
                                                                 simple_anno_size = unit(3, "mm"),
                                                                 col = c("chr3-5" = "black", "None" = "white")),
                                 Translocation.t32 = anno_simple(x = top_col_anno_df$Translocation.t32,
                                                                 simple_anno_size = unit(3, "mm"),
                                                                 col = c("chr3-2" = "black", "None" = "white")),
                                 Translocation.other = anno_text(x = Translocation.t3_other_text,
                                                                 which = "column", rot = 0, 
                                                                 location = 0.5, just = "center", height = unit(2, "mm"),
                                                                 gp = gpar(fill = "white", col = "black", fontsize = 5)))

# make row split ----------------------------------------------------------
rowsplit_vec <- ifelse(data_row_names.keep %in% c(vhl_complex_genes, "HIF1A", "EPAS1"), "VHL-HIF complex genes", "HIF targets")
rowsplit_factor <- factor(x = rowsplit_vec, levels = c("VHL-HIF complex genes", "HIF targets"))

# make row annotation -----------------------------------------------------
is_hif1a_target <- data_row_names.keep %in% hiftargets_df$target_genesymbol[hiftargets_df$source_genesymbol == "HIF1A"]
is_hif1a_target
colors_is_hif1a_target <- c("black", "white")
names(colors_is_hif1a_target) <- c("TRUE", "FALSE")

is_epas1_target <- data_row_names.keep %in% hiftargets_df$target_genesymbol[hiftargets_df$source_genesymbol == "EPAS1"]
is_epas1_target
colors_is_epas1_target <- c("black", "white")
names(colors_is_epas1_target) <- c("TRUE", "FALSE")
row_anno = rowAnnotation(Is_HIF1A_Target = anno_simple(x = as.character(is_hif1a_target), col = colors_is_hif1a_target),
                         Is_HIF2A_Target = anno_simple(x = as.character(is_epas1_target), col = colors_is_epas1_target))

# make heatmap body ------------------------------------------------------------
heatmapbody_color_fun <- colorRamp2(c(quantile(plot_data_mat, 0.025, na.rm=T), 
                                      quantile(plot_data_mat, 0.5, na.rm=T), 
                                      quantile(plot_data_mat, 0.975, na.rm=T)),
                                    c("blue", "white", "red"))

p <- Heatmap(matrix = plot_data_mat,
             col = heatmapbody_color_fun, 
             cluster_rows = T, row_split = rowsplit_factor, show_row_dend = F,
             row_title_rot = 0,
             cluster_columns = F, column_split = columnsplit_factor, show_column_dend = F,
             column_title_rot = 90, column_title_side = "bottom",
             right_annotation = row_anno, 
             top_annotation = top_col_anno,
             column_labels = ids_aliquot_wu, column_names_side = "top",
             show_heatmap_legend = F)
p
# make legend -------------------------------------------------------------
## make legend for top annotation
annotation_lgd = list(
  Legend(col_fun = heatmapbody_color_fun, 
         title = "Average expression value\nby cell type (Normalized)", 
         direction = "vertical"),
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
         title = "Bulk VHL Promoter Methylation"))

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "hifpathway.", "expression.", "tumor_vs_normal_epithelial_cells.", run_id, ".png")
png(filename = file2write, width = 2000, height = 2000, res = 150)
draw(object = p, 
     annotation_legend_side = "right", annotation_legend_list = annotation_lgd)
dev.off()
