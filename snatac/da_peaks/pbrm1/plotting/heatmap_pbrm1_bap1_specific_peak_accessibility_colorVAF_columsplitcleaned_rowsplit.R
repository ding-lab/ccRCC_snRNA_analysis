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
## input accessibility data
acc_pbrm1_down_df=fread(data.table = F, input = './Resources/snATAC_Processed_Data/Differential_Peaks/PBRM1_Specific/Accessibility/DOWN_PBRM1mutants_vs_nonMutans.Accessibility.20210623.tsv')
acc_pbrm1_up_df=fread(data.table = F, input ='./Resources/snATAC_Processed_Data/Differential_Peaks/PBRM1_Specific/Accessibility/UP_PBRM1mutants_vs_nonMutans.Accessibility.20210623.tsv')
acc_bap1_down_df=fread(data.table = F, input = './Resources/snATAC_Processed_Data/Differential_Peaks/BAP1_Specific/Accessibility/DOWN_BAP1mutants_vs_nonMutans.Accessibility.20210624.tsv')
acc_bap1_up_df=fread(data.table = F, input ='./Resources/snATAC_Processed_Data/Differential_Peaks/BAP1_Specific/Accessibility/UP_BAP1mutants_vs_nonMutans.Accessibility.20210624.tsv')

## input sample category
mut_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/mutation/annotate_cptac_sample_by_pbrm1_bap1_mutation/20210412.v1/PBRM1_BAP1_Mutation_Status_By_Case.20210412.v1.tsv")

# make data matrix --------------------------------------------------------
## add row names
row.names(acc_pbrm1_down_df)=paste(acc_pbrm1_down_df$Group.1,acc_pbrm1_down_df$Group.2,sep='_')
row.names(acc_pbrm1_up_df)=paste(acc_pbrm1_up_df$Group.1,acc_pbrm1_up_df$Group.2,sep='_')
row.names(acc_bap1_down_df)=paste(acc_bap1_down_df$Group.1,acc_bap1_down_df$Group.2,sep='_')
row.names(acc_bap1_up_df)=paste(acc_bap1_up_df$Group.1,acc_bap1_up_df$Group.2,sep='_')
## make row names consistent
rownames_plot <- rownames(acc_pbrm1_down_df)
## cbind
plotdata_df <- cbind(acc_pbrm1_down_df, 
                     acc_pbrm1_up_df[rownames_plot, grepl(pattern = "chr", x = colnames(acc_pbrm1_up_df))],
                     acc_bap1_down_df[rownames_plot, grepl(pattern = "chr", x = colnames(acc_bap1_down_df))],
                     acc_bap1_up_df[rownames_plot, grepl(pattern = "chr", x = colnames(acc_bap1_up_df))]
)
plotdata_df <- plotdata_df[,!duplicated(colnames(plotdata_df))]
# ## remove the BAP1-associated peaks from the PBRM1 associated peaks
# peaks_overlap <- c(colnames(acc_pbrm1_down_df), colnames(acc_pbrm1_up_df))
# peaks_overlap <- peaks_overlap[!(peaks_overlap %in% c(colnames(acc_bap1_down_df), colnames(acc_bap1_up_df)))]
# plotdata_df <- plotdata_df[, !(colnames(plotdata_df) %in% peaks_overlap)]
## change group names
plotdata_df$Group.1=gsub('_','-',plotdata_df$Group.1)
plotdata_df$Group.2=gsub('_','.',plotdata_df$Group.2)
## change column names
# colnames(plotdata_df)=gsub('\\.','\\-',colnames(plotdata_df))
## make ID column
plotdata_df$ID=paste(plotdata_df$Group.1,plotdata_df$Group.2,sep='_')
## change row names
rownames(plotdata_df)=plotdata_df$ID
## filter by row names
plotdata_df <- plotdata_df[!(rownames(plotdata_df) %in% c('C3L-00088-N_Tumor','C3N-01200-N_Tumor')),]
colnames(plotdata_df)[2:3]=c('Sample','Cell_type')
plotdata_df <- plotdata_df[plotdata_df$Cell_type %in% c('Tumor','PT'),]
plotdata_df2 <- plotdata_df[, colnames(plotdata_df)[grepl(pattern = "chr", x = colnames(plotdata_df))]]
plotdata_mat <- scale(x = plotdata_df2)
rownames(plotdata_mat) <- plotdata_df$Sample
## re-order rows
rownames_ordered <- c("C3L-00416-T2","C3L-00908-T1", 
                      "C3N-01200-T1", "C3L-01313-T1", "C3N-00317-T1","C3N-00437-T1", "C3L-01287-T1",
                      "C3N-01213-T1","C3L-00079-T1", "C3L-00610-T1", "C3N-00242-T1","C3L-00790-T1","C3L-00583-T1", "C3L-00004-T1", "C3N-00733-T1","C3L-01302-T1", 
                      "C3L-00448-T1",  "C3L-00917-T1", "C3L-00088-T1", "C3L-00088-T2","C3L-00010-T1", "C3L-00026-T1", "C3N-00495-T1","C3L-00096-T1", 
                      "C3N-01200-N", "C3L-00088-N")
plotdata_mat <- plotdata_mat[rownames_ordered,]

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## colors for gene set scores
summary(as.numeric(unlist(plotdata_mat)))
color_red <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[1]
color_blue <- RColorBrewer::brewer.pal(n = 5, name = "Set1")[2]
colors_heatmapbody = colorRamp2(c(-2, 0, 2), 
                                c(color_blue, "white", color_red))
colors_heatmapbody = colorRamp2(seq(-2, 2, 0.4), 
                                rev(brewer.pal(n = 11, name = "PiYG")))
colors_heatmapbody = colorRamp2(seq(-2, 2, 0.4), 
                                rev(brewer.pal(n = 11, name = "BrBG")))
colors_bap1_vaf <- colorRamp2(c(0, 0.5), c("white smoke", '#984EA3'))
colors_pbrm1_vaf <- colorRamp2(c(0, 0.5), c("white smoke", '#FF7F00'))
colors_peaktype <- c("Up-regulated" = "#E41A1C", "Down-regulated" = "#377EB8", "other" = "white smoke")

# make row annotation -----------------------------------------------------
## get BAP1 mutated cases
row_anno_df=plotdata_df %>%
  select(Sample, Cell_type) %>%
  mutate(Case = gsub(x = Sample, pattern = "\\-[A-Z][0-9]", replacement = ""))
row_anno_df$BAP1_mutation <- mapvalues(x = row_anno_df$Case, from = mut_df$Case, to = as.vector(mut_df$BAP1))
row_anno_df$BAP1_mutation[row_anno_df$BAP1_mutation == row_anno_df$Case] <- ""
row_anno_df$PBRM1_mutation <- mapvalues(x = row_anno_df$Case, from = mut_df$Case, to = as.vector(mut_df$PBRM1))
row_anno_df$PBRM1_mutation[row_anno_df$PBRM1_mutation == row_anno_df$Case] <- ""
row_anno_df$mutation_type <- mapvalues(x = row_anno_df$Case, from = mut_df$Case, to = as.vector(mut_df$mutation_category_sim))
row_anno_df$mutation_type[row_anno_df$mutation_type == row_anno_df$Case] <- "PT"
row_anno_df <- row_anno_df %>%
  mutate(BAP1_Mut_VAF = as.numeric(str_split_fixed(string = BAP1_mutation, pattern = "\\(|\\)", n = 3)[,2])) %>%
  mutate(PBRM1_Mut_VAF = as.numeric(str_split_fixed(string = PBRM1_mutation, pattern = "\\(|\\)", n = 3)[,2])) %>%
  mutate(BAP1_Mut_VAF = ifelse(is.na(BAP1_Mut_VAF), 0, BAP1_Mut_VAF)) %>%
  mutate(PBRM1_Mut_VAF = ifelse(is.na(PBRM1_Mut_VAF), 0, PBRM1_Mut_VAF))
rownames(row_anno_df) <- row_anno_df$Sample
row_anno_df <- row_anno_df[rownames_ordered,]
row_ha= rowAnnotation(#Cell_type=row_anno_df$Cell_type, 
  BAP1_Mut_VAF=row_anno_df$BAP1_Mut_VAF,
  PBRM1_Mut_VAF=row_anno_df$PBRM1_Mut_VAF,
  col=list(BAP1_Mut_VAF=colors_bap1_vaf,
           PBRM1_Mut_VAF=colors_pbrm1_vaf,
           Cell_type=c('PT'='#1B9E77','Tumor'='#E7298A')), 
  annotation_width = unit(2, "mm"), show_legend = F)

# make column annotation --------------------------------------------------
# column_ha <- HeatmapAnnotation(Is_BAP1_down_peak = anno_simple(x = as.character(colnames(plotdata_mat) %in% colnames(acc_bap1_down_df)), col = c("TRUE" = "purple", "FALSE" = "white smoke")),
#                                Is_BAP1_up_peak = anno_simple(x = as.character(colnames(plotdata_mat) %in% colnames(acc_bap1_up_df)), col = c("TRUE" = "purple", "FALSE" = "white smoke")),
#                                Is_PBRM1_down_peak = anno_simple(x = as.character(colnames(plotdata_mat) %in% colnames(acc_pbrm1_down_df)), col = c("TRUE" = "orange", "FALSE" = "white smoke")),
#                                Is_PBRM1_up_peak = anno_simple(x = as.character(colnames(plotdata_mat) %in% colnames(acc_pbrm1_up_df)), col = c("TRUE" = "orange", "FALSE" = "white smoke")), 
#                                annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 13))
"#377EB8"
"#E41A1C"
is_bap1_peaks_vec <- ifelse(colnames(plotdata_mat) %in% colnames(acc_bap1_down_df), "Down-regulated",
                                         ifelse(colnames(plotdata_mat) %in% colnames(acc_bap1_up_df), "Up-regulated", "other"))
is_pbrm1_peaks_vec <- ifelse(colnames(plotdata_mat) %in% colnames(acc_pbrm1_down_df), "Down-regulated",
                            ifelse(colnames(plotdata_mat) %in% colnames(acc_pbrm1_up_df), "Up-regulated", "other"))
column_ha <- HeatmapAnnotation(BAP1_associated_peak = anno_simple(x = is_bap1_peaks_vec, col = colors_peaktype[is_bap1_peaks_vec], height = unit(0.6, "cm")),
                               PBRM1_associated_peak = anno_simple(x = is_pbrm1_peaks_vec, col = colors_peaktype[is_pbrm1_peaks_vec],  height = unit(0.6, "cm")),
                               annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 17))

# make column split -------------------------------------------------------
column_split_vec <- ifelse(colnames(plotdata_mat) %in% colnames(acc_pbrm1_down_df), "PBRM1 down peak",
                           ifelse(colnames(plotdata_mat) %in% colnames(acc_pbrm1_up_df), "PBRM1 up peak", 
                                  ifelse(colnames(plotdata_mat) %in% colnames(acc_bap1_down_df), "BAP1 down peak", "BAP1 up peak")))

# make row split -------------------------------------------------------
row_split_vec <- row_anno_df$mutation_type
row_split_vec[row_split_vec == "Both mutated"] <- "BAP1&PBRM1\nmutated"
row_split_factor <- factor(row_split_vec, levels = c("BAP1&PBRM1\nmutated", "BAP1 mutated", 'PBRM1 mutated', "Non-mutants", 'PT'))

# plot --------------------------------------------------------------------
# p=ComplexHeatmap::Heatmap(matrix = plotdata_mat, col = colors_heatmapbody, name = "Relative peak\naccessibility", 
#                           ## rows
#                           show_row_names = T,  row_names_gp = gpar(fontsize = 13),
#                           show_row_dend=FALSE,  right_annotation=row_ha, row_split = row_split_factor, cluster_row_slices=F, 
#                           ## columns
#                           show_column_names = FALSE, show_column_dend=FALSE, column_title = NULL, top_annotation = column_ha,
#                           column_split = column_split_vec, cluster_column_slices=F,
#                           ## other
#                           show_heatmap_legend = F, use_raster = F)
list_lgd = list(
  Legend(title = "Peak accessibility direction vs. non-mutants", 
         title_gp = gpar(fontsize = 12),
         labels = names(colors_peaktype), legend_gp = gpar(fill = colors_peaktype), grid_width = unit(0.5, "cm"), grid_height = unit(0.75, "cm"),
         labels_gp =  gpar(fontsize = 12),direction = "horizontal", nrow = 1),
  Legend(col_fun = colors_heatmapbody, 
         title = "Relative peak\naccessibility",
         title_gp = gpar(fontsize = 12),
         labels_gp = gpar(fontsize = 12),
         legend_width = unit(3, "cm"), 
         direction = "horizontal"),
  Legend(col_fun = colors_bap1_vaf, 
         title = "BAP1 mutation VAF", 
         title_gp = gpar(fontsize = 12),
         labels_gp = gpar(fontsize = 12),
         legend_width = unit(3, "cm"),
         direction = "horizontal"),
  Legend(col_fun = colors_pbrm1_vaf, 
         title = "PBRM1 mutation VAF", 
         title_gp = gpar(fontsize = 12),
         labels_gp = gpar(fontsize = 12),
         legend_width = unit(3, "cm"),
         direction = "horizontal"))

# file2write <- paste0(dir_out, "PBRM1_specific_peak_accessibility.pdf")
# pdf(file2write,width=6.4, height=6.5, useDingbats = F)
# draw(object = p,
#      annotation_legend_side = "right", annotation_legend_list = list_lgd)
# dev.off()

p=ComplexHeatmap::Heatmap(matrix = plotdata_mat, col = colors_heatmapbody, name = "Peak\naccessibility", 
                          ## rows
                          show_row_names = F,  row_names_gp = gpar(fontsize = 13),
                          show_row_dend=FALSE,  left_annotation=row_ha, row_split = row_split_factor, row_title_rot = 0, row_title_gp = gpar(fontsize = 17),
                          cluster_row_slices=F, cluster_rows = F,
                          ## columns
                          show_column_names = FALSE, show_column_dend=FALSE, column_title = NULL, top_annotation = column_ha,
                          column_split = column_split_vec, cluster_column_slices=F,
                          ## other
                          show_heatmap_legend = F, use_raster = T)
file2write <- paste0(dir_out, "PBRM1_specific_peak_accessibility_raster.pdf")
pdf(file2write, width = 10, height=6.5, useDingbats = F)
draw(object = p,
     annotation_legend_side = "bottom", annotation_legend_list = list_lgd)
dev.off()

# p=ComplexHeatmap::Heatmap(matrix = plotdata_mat, col = colors_heatmapbody, name = "Peak\naccessibility", 
#                           ## rows
#                           show_row_names = T,  row_names_gp = gpar(fontsize = 13),
#                           show_row_dend=FALSE,  left_annotation=row_ha, row_split = row_split_factor, row_title_rot = 0,
#                           cluster_row_slices=F, cluster_rows = F,
#                           ## columns
#                           show_column_names = FALSE, show_column_dend=FALSE, column_title = NULL, top_annotation = column_ha,
#                           column_split = column_split_vec, cluster_column_slices=F,
#                           ## other
#                           show_heatmap_legend = F, use_raster = F)
# file2write <- paste0(dir_out, "PBRM1_specific_peak_accessibility.pdf")
# pdf(file2write,width=12, height=6.5, useDingbats = F)
# draw(object = p,
#      annotation_legend_side = "right", annotation_legend_list = list_lgd)
# dev.off()

