# Yige Wu @WashU Jul 2020

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
## input protein data
pho_df <- fread("./Resources/Bulk_Processed_Data/phosphoproteome/6_CPTAC3_CCRCC_Phospho_abundance_phosphosite_protNorm=2_CB.tsv", data.table = F)
## input bulk meta data for the entire set
bulk_meta_tab <- fread("~/Box/Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/cptac-metadata.csv")
## input the peaks to TF motifs
deg2tfgene_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/annotate_deg/annotate_tumor_vs_normal_degs_to_tf_by_snatac/20200916.v1/DEGs_with_DA_peaks.20200916.v1.tsv")

# get genes to plot -------------------------------------------------------
phosphosites <- data.frame(gene = c(rep("FLT1", 5), 
                                    rep("KDR", 8),
                                    "NRP1"),
                           phosphosite = c("Y1213", "Y794", "Y1313", "Y1169", "T680",
                                           "Y951", "Y996", "Y1054", "Y1059", "Y1175", "Y1214", "S1235", "S984",
                                           "Y297"))


# make ids to plot --------------------------------------------------------
# get the aliquot IDs for bulk corresponding to the snRNA aliquots --------
normal_bulk_aliquot_ids2plot <- bulk_meta_tab$Specimen.Label[bulk_meta_tab$Set.A == "yes" & bulk_meta_tab$Type == "Normal"]
normal_bulk_aliquot_ids2plot

case_ids2plot <- mapvalues(x = normal_bulk_aliquot_ids2plot, from = bulk_meta_tab$Specimen.Label, to = bulk_meta_tab$Case.ID)
case_ids2plot

tumor_bulk_aliquot_ids2plot <- mapvalues(x = case_ids2plot, from = bulk_meta_tab$Case.ID[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"], to = as.vector(bulk_meta_tab$Specimen.Label[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"]))
tumor_bulk_aliquot_ids2plot

# make heatmap body -------------------------------------------------------
## set the column names for the identifiers
header_col_names_keep <- c("Gene", "phosphosite", "Index", "ReferenceIntensity")
pho_tmp_df <- pho_df %>%
  filter(Gene == "KDR")
plot_data_df <- pho_df %>%
  mutate(phosphosite = str_split_fixed(string = Index, pattern = "_", n = 7)[,7]) %>%
  filter(phosphosite != "")  %>%
  select(header_col_names_keep, tumor_bulk_aliquot_ids2plot, normal_bulk_aliquot_ids2plot)
## filter the phosphosite
row_numbers_filter <- sapply(1:nrow(phosphosites), function(i, phosphosites_df, phosphosite_data_df) {
  gene_tmp <- phosphosites_df$gene[i]
  phosphosite_tmp <- phosphosites_df$phosphosite[i]
  row_numbers <- which(phosphosite_data_df$Gene == gene_tmp & grepl(x = phosphosite_data_df$phosphosite, pattern = phosphosite_tmp))
  return(row_numbers)
}, phosphosites_df = phosphosites, phosphosite_data_df = plot_data_df)
row_numbers_filter <- unlist(row_numbers_filter)
row_numbers_filter
plot_data_df <- plot_data_df[row_numbers_filter,]

## add a column by combining gene and phosphosite for filtering
plot_data_df <- plot_data_df %>%
  mutate(gene_phosphosite = paste0(Gene, "_", phosphosite))
## remove duplicates
plot_data_df <- plot_data_df[!duplicated(plot_data_df$gene_phosphosite),]
## transform data frame to matrix
plot_data_ma1 <- plot_data_df %>%
  select(tumor_bulk_aliquot_ids2plot, normal_bulk_aliquot_ids2plot)
rownames(plot_data_ma1) <- plot_data_df$gene_phosphosite
## take out the reference intensity
plot_data_mat <- as.matrix(plot_data_ma1) - as.vector(plot_data_df$ReferenceIntensity)

## get ids
ids_aliquot <- colnames(mat2plot)
sampletypes <- ifelse(ids_aliquot %in% tumor_bulk_aliquot_ids2plot, "Tumor", "Normal")

# make colors -------------------------------------------------------------
colors_sampletype <- c("Tumor" = "red", "Normal" = "green")
## make color function for heatmap body colors
color_blue <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[2]
color_red <- RColorBrewer::brewer.pal(n = 3, name = "Set1")[1]
heatmapbody_col_fun = colorRamp2(c(-2, 
                                   0, 
                                   2), 
                                 c(color_blue, "white", color_red))

# make column annotation --------------------------------------------------
ca = HeatmapAnnotation(Sample_Type = sampletypes,
                       col = list(Sample_Type = colors_sampletype), 
                       show_legend = F)


# plot heatmap ------------------------------------------------------------
p <- Heatmap(matrix = plot_data_mat, 
             col = heatmapbody_col_fun,
             ## row
             cluster_rows = F, show_row_dend = F,
             row_names_gp = gpar(fontsize = 12), row_names_side = "left",
             ## column
             column_split = sampletypes,
             cluster_columns = F, show_column_dend = F, 
             column_title_gp = gpar(fontsize = 18),
             show_column_names = F, 
             show_heatmap_legend = F)

annotation_lgd = list(
  Legend(col_fun = heatmapbody_col_fun, 
         title = "Protein abundance value\nlog2Intensity\n(Sample-Reference)", 
         legend_width = unit(6, "cm"),
         direction = "horizontal"))


## save heatmap to file
file2write <- paste0(dir_out, "Heatmap_by_sample.", ".pdf")
pdf(file2write, width = 20, height = 5)
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = annotation_lgd)
dev.off()
file2write <- paste0(dir_out, "Heatmap_by_sample.",".png")
png(file2write, width = 1200, height = 500, res = 150)
draw(object = p, 
     annotation_legend_side = "bottom", annotation_legend_list = annotation_lgd)
dev.off()


