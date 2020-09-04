# Yige Wu @WashU March 2020
## running on local
## for plotting average expression of known pathogenic pathway genes for each tumor subclusters (manually grouped)
## VHL-HIF pathway
## chromatin remodeling and epigenetic regulators
## MTOR pathway

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/aes.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the average expression calculated (RNA)
avg.exp.mat <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/averageexpression/averageexpression_tumor_cells_by_manual_subcluster/20200325.v1/averageexpression_tumor_cells_by_manual_subcluster.20200325.v1.tsv", data.table = F)
avg.exp.mat <- avg.exp.mat %>%
  rename(gene = V1)
## input ccRCC downstream genes
ccrcc_gene_alt_downstream_genes_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/dependencies/write_ccrcc_genetic_event_downstream_genes/20200302.v2/ccRCC_Genetic_Event_Downstream_Genes.20200302.v2.tsv", data.table = F)
## input the variable gene list
variable_genes_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/30_aliquot_integration/findvariablefeatures/findvariablefeatures_tumor_cells/20200310.v1/findvariablefeatures_tumor_cells.20200310.v1.tsv", data.table = F)
## input id meta data
id_metadata_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/make_meta_data/20191105.v1/meta_data.20191105.v1.tsv", data.table = F)
## input fraction of tumor subcluster cells with certain CNVs
cnv_state_count_aliquots <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/copy_number/summarize_cnv_fraction/cnv_fraction_in_malignant_nephron_epithelium_per_cluster/20200225.v1/Fraction_of_Malignant_Nephron_Epithelium_with_CNA_by_Gene.Per_Tumor_Subcluster.20200225.v1.tsv", data.table = F)

# make data matrix for heatmap body ---------------------------------------
## format the column names to only aliquot id 
data_col_names <- colnames(avg.exp.mat)[-1]
data_col_names.changed <- str_split_fixed(string = data_col_names, pattern = "\\.", n = 2)[,2]
## rename the data frame
colnames(avg.exp.mat) <- c("gene", data_col_names.changed)
## reformat data frame to matrix
plot_data_df <- avg.exp.mat
plot_data_mat <- as.matrix(plot_data_df[,-1])
plot_data_mat %>% head()
## add row names
rownames(plot_data_mat) <- plot_data_df$gene
plot_data_mat %>% head()
### get aliquot ids and case ids
tumorsubcluster_ids <- colnames(plot_data_mat)
aliquot_ids <- str_split_fixed(string = tumorsubcluster_ids, pattern = "_", n = 2)[,1]
case_ids <- mapvalues(x = aliquot_ids, from = id_metadata_df$Aliquot.snRNA, to = as.vector(id_metadata_df$Case))

# make top column annotation --------------------------------------------------
## make annotation data frame with copy number profile first
top_col_anno_df_long <- cnv_state_count_aliquots %>%
  filter(chr_region == "3p") %>%
  mutate(tumor_exp_subcluster_name = paste0(aliquot, "_EC", tumor_subcluster)) %>%
  filter(cna_state < 1) %>%
  group_by(tumor_exp_subcluster_name, gene_symbol) %>%
  summarise(Fraction = sum(Fraction))
top_col_anno_df_wide <- dcast(data = top_col_anno_df_long, formula = tumor_exp_subcluster_name ~ gene_symbol, value.var = "Fraction")
top_col_anno_df_wide <- merge(x = data.frame(tumor_exp_subcluster_name = colnames(plot_data_mat)), 
                              y = top_col_anno_df_wide, by = c("tumor_exp_subcluster_name"), all.x = T)
top_col_anno_df <- top_col_anno_df_wide[,-1]
rownames(top_col_anno_df) <- top_col_anno_df_wide$tumor_exp_subcluster_name
top_col_anno_df <- top_col_anno_df[colnames(plot_data_mat),]
top_col_anno_df[is.na(top_col_anno_df)] <- 0
# top_col_anno_df <- top_col_anno_df[colnames(plot_data_mat),]
## make annotation data frame with mutation type
### filter maf file by case id
maf_df <- maf_df %>%
  mutate(case_id = str_split_fixed(string = Tumor_Sample_Barcode, pattern = "_", n = 2)[,1]) %>%
  filter(case_id %in% case_ids)
var_class_df <- generate_somatic_mutation_matrix(pair_tab = ccRCC_SMGs, maf = maf_df)
var_class_anno_df <- t(var_class_df[,-1])
var_class_anno_df
var_class_anno_df <- var_class_anno_df[case_ids,]
var_class_anno_df
rownames(var_class_anno_df) <- rownames(plot_data_mat)


### make color for this annotation data frame
top_col_anno_colors <- lapply(colnames(top_col_anno_df), function(g) {
  if (g %in% c('3p', '5q', "14q")) {
    color_vector <- c("gain" = "red", "loss" = "blue")
  }
  if (g %in% ccRCC_SMGs) {
    color_vector <- variant_class_colors
  }
  return(color_vector)
})
names(top_col_anno_colors) <- colnames(top_col_anno_df)
### make top column annotation object
top_col_anno = HeatmapAnnotation(df = top_col_anno_df, show_legend = F)
top_col_anno = HeatmapAnnotation(Fraction_Cells_With_Deletion_by_gene = anno_barplot(x = top_col_anno_df, gp = gpar(fill = 1:4, col = 1:4)), show_legend = T)

# make bottom column annotation -------------------------------------------
bottom_col_anno = HeatmapAnnotation(foo = anno_text(case_ids, 
                                                    location = 0.5, just = "center",
                                                    gp = gpar(fill = uniq_case_colors[case_ids], col = "white", border = "black"),
                                                    width = max_text_width(case_ids)*1.2))

# plot all genes original values  ------------------------------------------------------
## make function for colors
heatmapbody_color_fun <- colorRamp2(c(quantile(plot_data_mat, 0.1, na.rm=T), 
                                      quantile(plot_data_mat, 0.5, na.rm=T), 
                                      quantile(plot_data_mat, 0.9, na.rm=T)),
                                    c("blue", "white", "red"))

p <- Heatmap(matrix = plot_data_mat,
             col = heatmapbody_color_fun,
             bottom_annotation = bottom_col_anno,
             # top_annotation = top_col_anno,
             show_heatmap_legend = T)

## make legend for top annotation
annotation_lgd = list(
  Legend(labels = c("BAP1", "VHL", "SETD2", "PBRM1"), 
         title = "Copy Number Deletion", 
         legend_gp = gpar(fill = 1:4, direction = "horizontal"),
         direction = "horizontal"))

## save heatmap
png(filename = paste0(dir_out, "avg_exp.tumorvariable_genes.by_tumor_manualsubcluster.heatmap.", run_id, ".png"), 
    width = 4500, height = 1500, res = 150)
### combine heatmap and heatmap legend
draw(object = p, 
     annotation_legend_side = "right", annotation_legend_list = annotation_lgd)
dev.off()

# plot all genes scaled across sample per gene   ------------------------------------------------------
## get the subset of genes
plot_genes_tmp <- as.vector(ccrcc_cna_genes_df$gene_symbol[ccrcc_cna_genes_df$chr_region == "3p"])
plot_genes_tmp <- intersect(plot_genes_tmp, avg.exp.mat$gene)
## get the subset of plot data
plot_data_mat_tmp <- plot_data_mat[plot_genes_tmp,]
## make function for colors
heatmapbody_color_fun <- colorRamp2(c(quantile(plot_data_mat_tmp, 0.1, na.rm=T), 
                                      quantile(plot_data_mat_tmp, 0.5, na.rm=T), 
                                      quantile(plot_data_mat_tmp, 0.9, na.rm=T)),
                                    c("blue", "white", "red"))
## make heatmap
p <- Heatmap(matrix = plot_data_mat_tmp,
             col = heatmapbody_color_fun,
             bottom_annotation = bottom_col_anno,
             top_annotation = top_col_anno,
             show_heatmap_legend = F)
p




