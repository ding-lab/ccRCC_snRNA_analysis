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
protein_tab <- fread("~/Box/Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/proteome/6_CPTAC3_CCRCC_Whole_abundance_gene_protNorm=2_CB.tsv", data.table = F)
## input bulk meta data for the entire set
bulk_meta_tab <- fread("~/Box/Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/cptac-metadata.csv")
## input filtered deg
genes_plot <- c("SLC16A9", "FTCD", "GLDC", "KCNJ15", "PDZD2", "RGN", "ABCC2")

# get the aliquot IDs for bulk corresponding to the snRNA aliquots --------
normal_bulk_aliquot_ids2plot <- bulk_meta_tab$Specimen.Label[bulk_meta_tab$Set.A == "yes" & bulk_meta_tab$Type == "Normal"]
normal_bulk_aliquot_ids2plot

case_ids2plot <- mapvalues(x = normal_bulk_aliquot_ids2plot, from = bulk_meta_tab$Specimen.Label, to = bulk_meta_tab$Case.ID)
case_ids2plot

tumor_bulk_aliquot_ids2plot <- mapvalues(x = case_ids2plot, from = bulk_meta_tab$Case.ID[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"], to = as.vector(bulk_meta_tab$Specimen.Label[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"]))
tumor_bulk_aliquot_ids2plot

# make heatmap body -------------------------------------------------------
## make the matrix to plot the heatmap
protein_tab2plot <- protein_tab %>%
  filter(Index %in% genes_plot) %>%
  select("Index", "ReferenceIntensity", tumor_bulk_aliquot_ids2plot, normal_bulk_aliquot_ids2plot)

protein_mat2plot <- protein_tab2plot %>%
  select(-Index) %>%
  select(-ReferenceIntensity)
rownames(protein_mat2plot) <- protein_tab2plot$Index

mat2plot <- as.matrix(protein_mat2plot) - as.vector(protein_tab2plot$ReferenceIntensity)
mat2plot <- mat2plot[genes_plot,]

# make column annotation --------------------------------------------------
ca = HeatmapAnnotation(Sample_Type = ifelse(colnames(mat2plot) %in% tumor_bulk_aliquot_ids2plot, "Tumor", "Normal"),
                       col = list(Sample_Type = c("Tumor" = "red", "Normal" = "green")))
# ra = rowAnnotation(Sample_Type = ifelse(rownames(mat2plot) %in% tumor_bulk_aliquot_ids2plot, "Tumor", "Normal"),
#                        col = list(Sample_Type = c("Tumor" = "red", "Normal" = "green")))

# plot heatmap ------------------------------------------------------------
p <- Heatmap(mat2plot,
             top_annotation = ca,
             cluster_columns = T, show_column_dend = F, 
             column_title = "Bulk Protein Level of Pseudotime Associated DEGs Between\nTumor and Normal Adjacent Tissue", column_title_gp = gpar(fontsize = 18),
             cluster_rows = F, row_names_gp = gpar(fontsize = 12), 
             name = "log2Intensity\n(Sample-Reference)", show_column_names = F)

## save heatmap to file
file2write <- paste0(dir_out, "Heatmap_by_sample.", ".pdf")
pdf(file2write, width = 20, height = 5)
print(p)
dev.off()
file2write <- paste0(dir_out, "Heatmap_by_sample.",".png")
png(file2write, width = 2000, height = 500, res = 150)
print(p)
dev.off()


