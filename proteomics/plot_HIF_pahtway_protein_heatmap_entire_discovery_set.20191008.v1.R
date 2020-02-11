# Yige Wu @ WashU 2019 Nov
## plot a heatmap with proteomics data from the discovery set data freeze

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")


# set run id --------------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)


# set output directory ----------------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input protein data ------------------------------------------------------
protein_tab <- fread("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/proteome/6_CPTAC3_CCRCC_Whole_abundance_gene_protNorm=2_CB.tsv", data.table = F)

# input bulk meta data ----------------------------------------------------
bulk_meta_tab <- fread("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/cptac-metadata.csv")

# specify the samples to be plotted ---------------------------------------
snRNA_aliquot_ids2plot <- c("CPT0019130004", "CPT0001260013", "CPT0086350004", "CPT0010110013", "CPT0001180011", "CPT0020120013", "CPT0001220012", "CPT0014450005")
  
# get the aliquot IDs for bulk corresponding to the snRNA aliquots --------
normal_bulk_aliquot_ids2plot <- bulk_meta_tab$Specimen.Label[bulk_meta_tab$Set.A == "yes" & bulk_meta_tab$Type == "Normal"]
normal_bulk_aliquot_ids2plot

case_ids2plot <- mapvalues(x = normal_bulk_aliquot_ids2plot, from = bulk_meta_tab$Specimen.Label, to = bulk_meta_tab$Case.ID)
case_ids2plot

tumor_bulk_aliquot_ids2plot <- mapvalues(x = case_ids2plot, from = bulk_meta_tab$Case.ID[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"], to = as.vector(bulk_meta_tab$Specimen.Label[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"]))
tumor_bulk_aliquot_ids2plot

# input TF table ----------------------------------------------------------
tf_tab <- fread(input = "./Ding_Lab/Projects_Current/TP53_shared_data/resources/PPI/TF_interactions.txt", data.table = F)
HIF_tf_tab <- tf_tab %>%
  filter(source_genesymbol %in% c("HIF1A", "EPAS1"))

# plot the HIF1alpha and HIF2alpha genes -------------------------------------
## specify the genes to be plotted 
genes2plot <- c("HIF1A", "EPAS1", "ARNT")

## make the matrix to plot the heatmap
protein_tab2plot <- protein_tab %>%
  filter(Index %in% genes2plot) %>%
  select("Index", "ReferenceIntensity", tumor_bulk_aliquot_ids2plot, normal_bulk_aliquot_ids2plot)

protein_mat2plot <- protein_tab2plot %>%
  select(-Index) %>%
  select(-ReferenceIntensity)
rownames(protein_mat2plot) <- protein_tab2plot$Index

mat2plot <- as.matrix(protein_mat2plot) - as.vector(protein_tab2plot$ReferenceIntensity)
mat2plot <- mat2plot[intersect(genes2plot, rownames(mat2plot)),]
## make column annotation
ca = HeatmapAnnotation(Sample_Type = ifelse(colnames(mat2plot) %in% tumor_bulk_aliquot_ids2plot, "Tumor", "Normal"),
                       col = list(Sample_Type = c("Tumor" = "red", "Normal" = "green")))

## plot heatmap
p <- Heatmap(mat2plot,
             column_names_gp = gpar(fontsize = 5),
             name = "log2Intensity\n(Sample-Reference)", 
             top_annotation = ca,
             cluster_columns = F,
             cluster_rows = F)

## save heatmap to file
file2write <- paste0(dir_out, "HIF_Protein_Expression.", run_id, ".png")
png(file2write, width = 2000, height = 400, res = 150)
print(p)
dev.off()


# plot the genes refined by snRNA-seq -------------------------------------
## specify the genes to be plotted 
genes2plot <- intersect(HIF_tf_tab$target_genesymbol, protein_tab$Index)
genes2plot <- c("VHL", "HIF1A", "EPAS1", genes2plot)

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

## make column annotation
ca = HeatmapAnnotation(Sample_Type = ifelse(colnames(mat2plot) %in% tumor_bulk_aliquot_ids2plot, "Tumor", "Normal"),
                       col = list(Sample_Type = c("Tumor" = "red", "Normal" = "green")))
# ra = rowAnnotation(Sample_Type = ifelse(rownames(mat2plot) %in% tumor_bulk_aliquot_ids2plot, "Tumor", "Normal"),
#                        col = list(Sample_Type = c("Tumor" = "red", "Normal" = "green")))

## plot heatmap
p <- Heatmap(mat2plot,
             name = "log2Intensity\n(Sample-Reference)", show_column_names = F,
             top_annotation = ca,
             cluster_columns = F,
             cluster_rows = F)

## save heatmap to file
file2write <- paste0(dir_out, "All_HIF_Downstream_Protein_Expression.", run_id, ".png")
png(file2write, width = 2000, height = 1000, res = 150)
print(p)
dev.off()


