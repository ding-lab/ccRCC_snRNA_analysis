# Yige Wu @WashU March 2020
## for each individual sample tumor cell reclustered, plot UMAP

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
## input protein data
rna_df <- fread("./Resources/Bulk_Processed_Data/mRNA/RNA_rpkm_tumor_normal.tsv", data.table = F)
## input bulk meta data for the entire set
bulk_meta_df <- fread(data.table = F, input = "./Resources/Bulk_Processed_Data/Meta_Data/cptac-metadata.csv")
## input mutation classificatin
mut_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/mutation/annotate_cptac_sample_by_pbrm1_bap1_mutation/20210310.v1/PBRM1_BAP1_Mutation_Status_By_Case.20210310.v1.tsv")

# specify the genes to be plotted  -------------------------------------
genes2plot <- c("EPHA6", "ABCC3", "FTO", "COL23A1", "CA9", "SEMA6A", "NDRG1", "EGFR", "CP")
bulk_meta_df <- bulk_meta_df %>%
  mutate(easyid = ifelse(!is.na(Type), ifelse(Type == "Tumor", paste0(Case.ID, "-T1"), paste0(Case.ID, "-N")), Case.ID))

# make heatmap body -------------------------------------------------------
## filter
rna_filtered_df <- rna_df %>%
  filter(geneID %in% genes2plot)
## transform
rna_filtered_mat <- as.matrix(rna_filtered_df[, -1])
plot_raw_data_mat <- log2(rna_filtered_mat+1)
## filter
plot_raw_data_mat <- plot_raw_data_mat[, bulk_meta_df$RNA.ID[!is.na(bulk_meta_df$Type) &  bulk_meta_df$Type == "Tumor"]]
## get ids
geneids_plot <- rna_filtered_df$geneID
rnaids_plot <- colnames(plot_raw_data_mat)
caseids_plot <- mapvalues(x = rnaids_plot, from = bulk_meta_df$RNA.ID, to = as.vector(bulk_meta_df$Case.ID))
easyids_plot <- mapvalues(x = rnaids_plot, from = bulk_meta_df$RNA.ID, to = as.vector(bulk_meta_df$easyid))
sampletypes_plot <- mapvalues(x = rnaids_plot, from = bulk_meta_df$RNA.ID, to = as.vector(bulk_meta_df$Type))
## scale by row
plot_data_mat <- t(apply(plot_raw_data_mat, 1, scale))
## add row names and colnames
rownames(plot_data_mat) <- geneids_plot
colnames(plot_data_mat) <- easyids_plot

# specify colors ----------------------------------------------------------
## specify color for NA values
color_na <- "grey50"
## make color function for heatmap body colors
color_blue <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[2]
color_red <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
summary(as.vector(plot_data_mat))
colors_heatmapbody = circlize::colorRamp2(c(-2,
                                            0,
                                            2),
                                          c(color_blue, "white", color_red))
## make colors for the original unscaled expression
# summary(as.vector(orig_avgexp_vec))
colors_unscaledexp = circlize::colorRamp2(1:9, 
                                          RColorBrewer::brewer.pal(name = "RdPu", n = 9))

# make column split -------------------------------------------------------
col_split_vec <- mapvalues(x = caseids_plot, from = mut_df$Case, to = as.vector(mut_df$mutation_category_sim))
col_split_vec[col_split_vec %in% caseids_plot] <- "Non-mutants"
col_split_vec[sampletypes_plot == "Normal"] <- "NAT"

# make row annotatin ------------------------------------------------------
orig_avgexp_vec <- rowMeans(x = plot_raw_data_mat, na.rm = T)
row_anno_obj <- rowAnnotation(Unscaled_Expression = anno_simple(x = orig_avgexp_vec, col = colors_unscaledexp), 
                              annotation_name_side = "top")

# plot heatmap ------------------------------------------------------------
p <- Heatmap(matrix = plot_data_mat, col = colors_heatmapbody,
             ## row
             # row_split = row_split_vec, show_row_dend = F,
             row_names_gp = gpar(fontsize = 10), row_title_rot = 0,
             right_annotation = row_anno_obj, 
             ## column
             column_km = 3, column_km_repeats = 100,
             # column_split = col_split_vec, show_column_dend = F, 
             # cluster_column_slices = T, column_title_rot = 90,
             show_column_names = F)

## save heatmap to file
file2write <- paste0(dir_out, "DEGs",  ".png")
png(file2write, width = 2000, height = 2000, res = 150)
print(p)
dev.off()
file2write <- paste0(dir_out, "DEGs",  ".pdf")
pdf(file2write, width = 15, height = 6)
print(p)
dev.off()


