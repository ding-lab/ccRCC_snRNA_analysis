# Yige Wu @ WashU 2020 Feb
## plot a heatmap with CNV status for the frequently altered genes and their expression

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input RNA expression
rna_tab <- fread("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/mRNA/RNA_rpkm_tumor_normal.tsv", data.table = F)
## input CNA matrix
cna_tab <- loadCNAstatus()
## input bulk meta data
bulk_meta_tab <- fread("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/cptac-metadata.csv")
## input snRNA sample set
# snRNA_meta_data <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/make_meta_data/20191105.v1/meta_data.20191105.v1.tsv", data.table = F)
srat_paths <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/individual_cluster/write_path_to_seurat_objects_on_box/20200219.v1/Seurat_Object_Paths.20200219.v1.tsv", data.table = F)
## set NonNA threshold
num_nonna <- 0

# set the aliquot IDs for bulk corresponding to the snRNA aliquots --------
case_ids2plot <- unique(srat_paths$Case)
case_ids2plot

normal_bulk_aliquot_ids2plot <- mapvalues(x = case_ids2plot, from = bulk_meta_tab$Case.ID[bulk_meta_tab$Type == "Normal" & bulk_meta_tab$Set.A == "yes"], to = as.vector(bulk_meta_tab$RNA.ID[bulk_meta_tab$Type == "Normal" & bulk_meta_tab$Set.A == "yes"]))
normal_bulk_aliquot_ids2plot <- normal_bulk_aliquot_ids2plot[!(normal_bulk_aliquot_ids2plot %in% case_ids2plot)]
normal_bulk_aliquot_ids2plot

tumor_bulk_aliquot_ids2plot <- mapvalues(x = case_ids2plot, from = bulk_meta_tab$Case.ID[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"], to = as.vector(bulk_meta_tab$RNA.ID[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"]))
tumor_bulk_aliquot_ids2plot

# plot gene expression data -------------------------------------
## specify the genes to be plotted 
genes2plot <- ccrcc_cna_genes_df$gene_symbol[ccrcc_cna_genes_df$chr_region %in% c("3p", "5q", "14q")]
## filter the CNVs
cna_tab <- cna_tab[cna_tab$gene %in% genes2plot, c("gene", unique(srat_paths$Case))]
## make the matrix to plot the heatmap
exp_tab2plot <- rna_tab %>%
  filter(geneID %in% genes2plot) %>%
  select("geneID", tumor_bulk_aliquot_ids2plot, normal_bulk_aliquot_ids2plot)

exp_mat2plot <- exp_tab2plot %>%
  select(-geneID)
exp_mat2plot <- as.matrix(exp_mat2plot)
rownames(exp_mat2plot) <- exp_tab2plot$geneID
mat2plot <- log10(exp_mat2plot)
mat2plot <- mat2plot[intersect(genes2plot, rownames(mat2plot)),]
colnames(mat2plot) <- c(paste0(mapvalues(x = tumor_bulk_aliquot_ids2plot, from = bulk_meta_tab$RNA.ID, to = bulk_meta_tab$Case.ID), "_Tumor"), 
                        paste0(mapvalues(x = normal_bulk_aliquot_ids2plot, from = bulk_meta_tab$RNA.ID, to = bulk_meta_tab$Case.ID), "_Normal"))
## make column annotation
rownames(cna_tab) <- cna_tab$gene
col_anno_df <- t(cna_tab[,-1])
col_anno_df <- col_anno_df[,as.vector(genes2plot)]
col_anno_df <- rbind(col_anno_df, matrix(data = NA, nrow = length(normal_bulk_aliquot_ids2plot), ncol = ncol(col_anno_df)))
col_anno_df <- as.data.frame(col_anno_df)
### make color annotation
col_anno_colors <- lapply(genes2plot, function(g) {
  color_vector <- c("amplification" = "red", "deletion" = "blue", "neutral" = "white")
  return(color_vector)
})
names(col_anno_colors) <- genes2plot
bca = HeatmapAnnotation(df = col_anno_df, col = col_anno_colors)
tca = HeatmapAnnotation(Sample_Type = ifelse(grepl(x = colnames(mat2plot), pattern = "Tumor"), "Tumor", "Normal"),
                        col = list(Sample_Type = c("Tumor" = "red", "Normal" = "green")))
## get color corresponding to values
col_rna <- colorRamp2(c(quantile(mat2plot, 0.1, na.rm=T), mean(mat2plot, na.rm=T), quantile(mat2plot, 0.9, na.rm=T)),c("blue", "white", "red"))
# col_rna <- colorRamp2(c(-1.5, 0, 1.5),c("blue", "white", "red"))

## plot heatmaps
p <- Heatmap(mat2plot[(rowSums(!is.na(mat2plot)) >= num_nonna),],
             col = col_rna, na_col = "grey",
             name = "log10(RPKM)", 
             row_names_gp = gpar(fontsize = 10),
             column_names_gp = gpar(fontsize = 10),
             bottom_annotation = bca,
             top_annotation = tca,
             cluster_columns = F,
             cluster_rows = F)
p

## save heatmap to PNG
file2write <- paste0(dir_out, "CNA_Genes_CNV_and_mRNA_Expression.", run_id, ".png")
png(file2write, width = 1400, height = 700, res = 150)
print(p)
dev.off()

## save heatmap to PDF
file2write <- paste0(dir_out, "CNA_Genes_CNV_and_mRNA_Expression.", run_id, ".pdf")
pdf(file2write, width = 12, height = 5)
print(p)
dev.off()
