# Yige Wu @WashU March 2020
## running on local

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/load_data.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## load bulk protein data
protein_df <- loadParseProteomicsData(expression_type = "protein")
## input bulk meta data
metadata_bulk_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Meta_Data/cptac-metadata.csv")
## input id meta data
metadata_sn_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210305.v1/meta_data.20210305.v1.tsv", data.table = F)
## input mutation file
maf_df <- loadMaf()
## input te bulk genomics/methylation events
bulk_sn_omicsprofile_df <- fread(input = "./Resources/Analysis_Results/data_summary/merge_bulk_sn_profiles/20200512.v1/bulk_sn_omics_profile.20200512.v1.tsv", data.table = F)

# preprocess --------------------------------------------------------------
## set parameters
exp_genes <- unique(c("PBRM1", "BAP1","SETD2", 
                      "KDM5C", "KDM6A",
                      pbaf_genes))
## preprocess data
metadata_bulk_df <- metadata_bulk_df %>%
  mutate(easy_id = ifelse(Type == "Tumor", paste0(Case.ID, "-T1"), paste0(Case.ID, "-N")))
### preprocess the protein data
protein_df <- protein_df %>%
  filter(Index %in% exp_genes)
rownames(protein_df) <- protein_df$Index
protein_df <- protein_df[exp_genes[exp_genes %in% protein_df$Index],]
### merge meta data
metadata_merged_df <- merge(x = metadata_bulk_df %>%
                              filter(Type == "Tumor"), 
                            y = metadata_sn_df %>%
                              filter(Sample_Type == "Tumor"), 
                            by.x = c("easy_id", "Case.ID"), by.y = c("Aliquot.snRNA.WU", "Case"), all = T)
easyids_plot <- metadata_merged_df$easy_id
### make mutation matrix
mut_mat <- get_mutation_class_sim_matrix(pair_tab = ccRCC_drivers, maf = maf_df)
mut_df <- t(mut_mat[,-1]) %>% as.data.frame()
mut_df$Case <- rownames(mut_df)

# merge plot data --------------------------------------------------------------
## scale protein data
colnames_protein <- colnames(protein_df)
aliquotids_protein <- colnames_protein[grepl(pattern = "CPT", x = colnames_protein)]
protein_scaled_df <- protein_df[,aliquotids_protein]
protein_scaled_df <- protein_scaled_df - protein_df$ReferenceIntensity
protein_scaled_mat <- as.matrix(protein_scaled_df)
rownames(protein_scaled_mat) <- protein_df$Index
### rename column names
easyids_protein <- mapvalues(x = aliquotids_protein, from = metadata_bulk_df$Specimen.Label, to = as.vector(metadata_bulk_df$easy_id))
colnames(protein_scaled_mat) <- easyids_protein
# easyids_plot <- easyids_protein[grepl(pattern = "T1", x = easyids_protein)]
## add NA data for T2 and T3 segments
easyids_pro_na <- easyids_plot[!(easyids_plot %in% easyids_protein)]
protein_na_mat <- matrix(data = 0, nrow = nrow(protein_scaled_mat), ncol = length(easyids_pro_na))
colnames(protein_na_mat) <- easyids_pro_na
rownames(protein_na_mat) <- rownames(protein_scaled_mat)
protein_scaled_padded_mat <- cbind(protein_na_mat, protein_scaled_mat)
protein_scaled_filtered_mat <- protein_scaled_padded_mat[,easyids_plot]
## scale RNA data

## map ids
## merge protein and RNA data
plot_data_mat <- protein_scaled_filtered_mat

# make colors -------------------------------------------------------------
## make colors for the heatmap body
color_blue <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[2]
color_red <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
summary(as.vector(plot_data_mat))
colors_heatmapbody = colorRamp2(c(-0.5, 
                                  0, 
                                  0.5), 
                                c(color_blue, "white", color_red))
color_na <- "grey70"
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
top_col_anno_df <- merge(x = metadata_merged_df,
                         y = mut_df, 
                         by.x = c("Case.ID"), by.y = c("Case"), all.x = T)
rownames(top_col_anno_df) <- top_col_anno_df$easy_id
top_col_anno_df <- top_col_anno_df[easyids_plot,]
top_col_anno_df <- top_col_anno_df %>%
  mutate(Sample_atWashU = as.character(!is.na(Aliquot.snRNA))) %>%
  mutate(snRNA_available = as.character(!is.na(snRNA_available) & snRNA_available == TRUE)) %>%
  mutate(snATAC_available = as.character(!is.na(snATAC_available) & snATAC_available == TRUE))
### make text for translocations
top_col_anno = HeatmapAnnotation(
  snRNA_Availability = anno_simple(x = top_col_anno_df$snRNA_available,
                                   simple_anno_size = unit(5, "mm"),
                                   gp = gpar(col = color_gridline), col = colors_sn_data),
  snATAC_Availability = anno_simple(x = top_col_anno_df$snATAC_available,
                                    simple_anno_size = unit(5, "mm"), 
                                    gp = gpar(col = color_gridline),  col = colors_sn_data),
  Sample_atWashU = anno_simple(x = top_col_anno_df$Sample_atWashU,
                               gp = gpar(col = color_gridline),
                               simple_anno_size = unit(5, "mm"),
                               col = c("TRUE" = "black", "FALSE" = "grey80")),
  PBRM1_Mutation = anno_simple(x = top_col_anno_df$PBRM1,
                               gp = gpar(col = color_gridline),
                               simple_anno_size = unit(5, "mm"),
                               col = colors_variant_class_sim),
  BAP1_Mutation = anno_simple(x = top_col_anno_df$BAP1,
                              gp = gpar(col = color_gridline),
                              simple_anno_size = unit(5, "mm"),
                              col = colors_variant_class_sim),
  VHL_Mutation = anno_simple(x = top_col_anno_df$VHL,
                             gp = gpar(col = color_gridline), 
                             simple_anno_size = unit(5, "mm"),
                             col = colors_variant_class_sim),
  annotation_name_side = "left", annotation_name_gp = gpar(fontsize = 10))


# plot heatmap body with white-yellow-red ------------------------------------------------------
## make heatmap
p <- Heatmap(matrix = plot_data_mat, col = colors_heatmapbody,
             na_col = color_na,
             ## column
             cluster_columns = T,
             # column_split = col_split_factor,
             # column_order = order(match(plot_data_df$Aliquot_Suffix, c("T1", "T2", "T3", "N"))),
             # bottom_annotation = bottom_col_anno, 
             top_annotation = top_col_anno,
             column_gap = unit(x = 0, units = "mm"), show_column_dend = F, column_names_gp = gpar(fontsize = 10), column_names_side = "top",
             ## row
             show_row_dend = F, row_names_side = "left", cluster_rows = F,
             show_heatmap_legend = F)
p
## make legend for top annotation
annotation_lgd = list(
  Legend(labels = names(colors_variant_class_sim), labels_gp = gpar(fontsize = 9),
         title = "Somatic Mutation Status", title_gp = gpar(fontsize = 9),
         legend_gp = gpar(fill = colors_variant_class_sim), border = color_gridline),
  Legend(col_fun = colors_heatmapbody, 
         title = "Protein abundance value\nlog2Intensity\n(Sample-Reference)"))
# ## save heatmap as png
# png(filename = paste0(dir_out, "data_availability",".png"), 
#     width = 1600, height = 1600, res = 150)
# ### combine heatmap and heatmap legend
# draw(object = p, 
#      annotation_legend_side = "right", annotation_legend_list = annotation_lgd)
# dev.off()
## save heatmap as pdf
pdf(paste0(dir_out, "data_availability", ".pdf"), 
    width = 20, height = 6)
### combine heatmap and heatmap legend
draw(object = p, 
     annotation_legend_side = "right", annotation_legend_list = annotation_lgd)
dev.off()

