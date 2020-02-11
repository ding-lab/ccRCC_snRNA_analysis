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

# input meta data to distinguish tumor vs normal --------------------------
meta_tab <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/make_meta_data/20191105.v1/meta_data.20191105.v1.tsv", data.table = F)

# specify the samples to be plotted ---------------------------------------
snRNA_aliquot_ids2plot <- c("CPT0019130004", "CPT0001260013", "CPT0086350004", "CPT0010110013", "CPT0001180011", "CPT0020120013", "CPT0001220012", "CPT0014450005")
  
# get the aliquot IDs for bulk corresponding to the snRNA aliquots --------
case_ids2plot <- mapvalues(x = snRNA_aliquot_ids2plot, from = meta_tab$Aliquot.snRNA, to = as.vector(meta_tab$Case))
case_ids2plot

tumor_bulk_aliquot_ids2plot <- mapvalues(x = case_ids2plot, from = meta_tab$Case[meta_tab$Sample_Type == "Tumor" & meta_tab$Is_discovery_set == T], to = as.vector(meta_tab$Aliquot.bulk[meta_tab$Sample_Type == "Tumor" & meta_tab$Is_discovery_set == T]))
tumor_bulk_aliquot_ids2plot

normal_bulk_aliquot_ids2plot <- mapvalues(x = case_ids2plot, from = meta_tab$Case[meta_tab$Sample_Type == "Normal" & meta_tab$Is_discovery_set == T], to = as.vector(meta_tab$Aliquot.bulk[meta_tab$Sample_Type == "Normal" & meta_tab$Is_discovery_set == T]))
normal_bulk_aliquot_ids2plot

# specify the genes to be plotted -----------------------------------------
genes2plot <- c("HIF1A", "EPAS1", "MET", "CA9", "PGK1", "PKM", "LDHA", "PFK1", "VEGFA", "EGLN1", "EGLN3")

# make the matrix to plot the heatmap -------------------------------------
protein_tab2plot <- protein_tab %>%
  filter(Index %in% genes2plot) %>%
  select("Index", "ReferenceIntensity", tumor_bulk_aliquot_ids2plot, normal_bulk_aliquot_ids2plot)

protein_mat2plot <- protein_tab2plot %>%
  select(-Index) %>%
  select(-ReferenceIntensity)
rownames(protein_mat2plot) <- protein_tab2plot$Index

mat2plot <- as.matrix(protein_mat2plot) - as.vector(protein_tab2plot$ReferenceIntensity)
# protein_mat2plot <- t(protein_mat2plot)


# make column annotation --------------------------------------------------
ca = HeatmapAnnotation(Sample_Type = ifelse(colnames(mat2plot) %in% tumor_bulk_aliquot_ids2plot, "Tumor", "Normal"),
                       col = list(Sample_Type = c("Tumor" = "red", "Normal" = "green")))




# plot heatmap ------------------------------------------------------------
p <- Heatmap(mat2plot,
             name = "log2Intensity - Reference", 
             top_annotation = ca,
             # col = col_fun,
             cluster_columns = F,
             cluster_rows = F)
file2write <- paste0(dir_out, "HIF_Downstream_Protein_Expression.", run_id, ".png")
png(file2write, width = 1500, height = 500, res = 150)
print(p)
dev.off()


