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
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# input dependencies ------------------------------------------------------
## input protein data
protein_tab <- fread("~/Box/Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/proteome/6_CPTAC3_CCRCC_Whole_abundance_gene_protNorm=2_CB.tsv", data.table = F)
## input bulk meta data for the entire set
bulk_meta_tab <- fread("~/Box/Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/cptac-metadata.csv")
## input hif targets
hif_targets_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/write_hif_targets/20200428.v1/HIF_Target_Genes.20200428.v1.tsv")
## input DEG for each cell group
deg_df <- fread(input = "./Resources/Analysis_Results/findmarkers/findmarkers_by_celltype/findallmarker_wilcox_cellgroup_on_katmai/20200714.v2/findallmarkers_wilcox_bycellgroup.pos..logfcthreshold0.1.minpct0.1.mindiffpct0.1.tsv", data.table = F)

# get the aliquot IDs for bulk corresponding to the snRNA aliquots --------
normal_bulk_aliquot_ids2plot <- bulk_meta_tab$Specimen.Label[bulk_meta_tab$Set.A == "yes" & bulk_meta_tab$Type == "Normal"]
normal_bulk_aliquot_ids2plot

case_ids2plot <- mapvalues(x = normal_bulk_aliquot_ids2plot, from = bulk_meta_tab$Specimen.Label, to = bulk_meta_tab$Case.ID)
case_ids2plot

tumor_bulk_aliquot_ids2plot <- mapvalues(x = case_ids2plot, from = bulk_meta_tab$Case.ID[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"], to = as.vector(bulk_meta_tab$Specimen.Label[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"]))
tumor_bulk_aliquot_ids2plot

# input TF table ----------------------------------------------------------
tf_tab <- fread(input = "./Ding_Lab/Projects_Current/TP53_shared_data/resources/PPI/TF_interactions.txt", data.table = F)
hif_targets_df <- tf_tab %>%
  filter(source_genesymbol %in% c("HIF1A", "EPAS1"))


# specify the genes to be plotted  -------------------------------------
genes2plot <- intersect(hif_targets_df$target_genesymbol, protein_tab$Index)
genes2plot <- intersect(genes2plot, unique(deg_df$gene))
genes2plot <- unique(c("VHL", "HIF1A", "EPAS1", genes2plot))

# make heatmap body -------------------------------------------------------
## make the matrix to plot the heatmap
protein_tab2plot <- protein_tab %>%
  filter(Index %in% genes2plot) %>%
  select("Index", "ReferenceIntensity", tumor_bulk_aliquot_ids2plot, normal_bulk_aliquot_ids2plot)

protein_mat2plot <- protein_tab2plot %>%
  select(-Index) %>%
  select(-ReferenceIntensity)
rownames(protein_mat2plot) <- protein_tab2plot$Index

mat2plot <- as.matrix(protein_mat2plot) - as.vector(protein_tab2plot$ReferenceIntensity)
mat2plot <- mat2plot[genes2plot,]


# make column annotation --------------------------------------------------
ca = HeatmapAnnotation(Sample_Type = ifelse(colnames(mat2plot) %in% tumor_bulk_aliquot_ids2plot, "Tumor", "Normal"),
                       col = list(Sample_Type = c("Tumor" = "red", "Normal" = "green")))
# ra = rowAnnotation(Sample_Type = ifelse(rownames(mat2plot) %in% tumor_bulk_aliquot_ids2plot, "Tumor", "Normal"),
#                        col = list(Sample_Type = c("Tumor" = "red", "Normal" = "green")))

# make row split ----------------------------------------------------------
genes_multicelltypeexpr <- names(table(deg_df$gene)[table(deg_df$gene) == 2])
genes_tumorcellexpr <- deg_df$gene[deg_df$cluster == "Tumor cells" & !(deg_df$gene %in% genes_multicelltypeexpr)]
genes_normalepitheliumexpr <- deg_df$gene[deg_df$cluster == "Normal epithelial cells" & !(deg_df$gene %in% genes_multicelltypeexpr)]
genes_stromaexpr <- deg_df$gene[deg_df$cluster == "Stroma" & !(deg_df$gene %in% genes_multicelltypeexpr)]
genes_immuneexpr <- deg_df$gene[deg_df$cluster == "Immune" & !(deg_df$gene %in% genes_multicelltypeexpr)]
genes_other <- genes2plot[!(genes2plot %in% c(genes_tumorcellexpr, genes_normalepitheliumexpr, genes_stromaexpr, genes_immuneexpr))]
gene_celltype_exp_cat_df <- data.frame(gene = c(genes_other, 
                                                genes_tumorcellexpr, 
                                                genes_normalepitheliumexpr,
                                                genes_stromaexpr, 
                                                genes_immuneexpr),
                                       gene_celltypeexp_cat = c(rep("Other", length(genes_other)),
                                                                rep("TumorCells\nExpressed", length(genes_tumorcellexpr)),
                                                                rep("NormalEpithelium\nExpressed", length(genes_normalepitheliumexpr)),
                                                                rep("Stroma\nExpressed", length(genes_stromaexpr)),
                                                                rep("Immune\nExpressed", length(genes_immuneexpr))))
row_split_vec <- mapvalues(x = rownames(mat2plot), from = gene_celltype_exp_cat_df$gene, to = as.vector(gene_celltype_exp_cat_df$gene_celltypeexp_cat))
row_split_factor <- factor(x = row_split_vec, levels = c("Other", "TumorCells\nExpressed", "NormalEpithelium\nExpressed", "Stroma\nExpressed", "Immune\nExpressed"))
# plot heatmap ------------------------------------------------------------
p <- Heatmap(mat2plot,
             top_annotation = ca,
             cluster_columns = T, show_column_dend = F,
             cluster_rows = F, row_names_gp = gpar(fontsize = 12), row_split = row_split_factor, row_title_gp = gpar(fontsize = 12),
             name = "log2Intensity\n(Sample-Reference)", show_column_names = F)

## save heatmap to file
file2write <- paste0(dir_out, "All_HIF_Downstream_Protein_Expression.", run_id, ".pdf")
pdf(file2write, width = 20, height = 7)
print(p)
dev.off()
# file2write <- paste0(dir_out, "All_HIF_Downstream_Protein_Expression.", run_id, ".png")
# png(file2write, width = 2000, height = 1000, res = 150)
# print(p)
# dev.off()


