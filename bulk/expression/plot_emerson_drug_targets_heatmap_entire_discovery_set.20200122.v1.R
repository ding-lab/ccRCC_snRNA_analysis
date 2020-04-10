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


# set NonNA threshold -----------------------------------------------------
num_nonna <- 40
row_fontsize <- 9

# input gene list ---------------------------------------------------------
apoptosis_genes <- read_excel("./Ding_Lab/Projects_Current/RCC/ccRCC_Drug/Resources/Drug_Lists/Kidney_GBM_Genes_Jan2020.xlsx", 
                              sheet = "Apoptosis")
antiapop_genes <- apoptosis_genes$Gene[apoptosis_genes$Role %in% c("Anti-apoptosis", "IAP")]
proapop_genes <- apoptosis_genes$Gene[apoptosis_genes$Role %in% c("Pro-apoptosis")]
## targets by the IAP inhibitors: https://www.medchemexpress.com/Targets/IAP.html?src=google-product&gclid=EAIaIQobChMIr9eU1IWY5wIVEPDACh1IGA-_EAAYAiAAEgKxLPD_BwE
# IAP (Inhibitors of Apoptosis) is a family of functionally and structurally related proteins, which serve as endogenous inhibitors of programmed cell death (apoptosis). A common feature of all IAPs is the presence of a BIR in one to three copies. The human IAP family consists of 8 members, and IAP homologs have been identified in numerous organisms. The members of the IAPs included IAPs, Cp-IAP, Op-IAP, XIAP, c-IAPl, C-IAP2, NAIP, Livin and Survivin. The best characterized IAP is XIAP, which binds caspase-9, caspase-3 and caspase 7, thereby inhibiting their activation and preventing apoptosis. Also cIAP1 and cIAP2 have been shown to bind caspases, although how the IAPs inhibit apoptosis mechanistically at the molecular level is not completely understood.
IAP_genes <- c("BIRC2", "BIRC3", "XIAP", "BIRC7", "BIRC5")
cdk_genes <- c("CDK4", "CDK6", "CDK7", "WEE1")
mtor_genes <- c("MTOR", "AKTS1", "DEPTOR", "RPTOR", "MLST8", "MAPKAP1", "PRR5", "RICTOR")
## targets by carbozantinib: https://www.drugbank.ca/drugs/DB08875
rtk_genes <- c("MET", "AXL", "FLT1", "KDR", "FLT3", "KIT", "RET", "NTRK2", "TEK")


# plot the HIF1alpha and HIF2alpha genes -------------------------------------
## specify the genes to be plotted 
genes2plot <- c(IAP_genes, cdk_genes, mtor_genes, rtk_genes)

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

## get color corresponding to values
col_protein = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))


## plot heatmaps
p_rtk <- Heatmap(mat2plot[(rownames(mat2plot) %in% rtk_genes) & (rowSums(!is.na(mat2plot)) >= num_nonna),],
                      row_title = "RTKs",
                      row_title_gp = gpar(fontsize = 12),
                      col = col_protein,
                      name = "log2Intensity\n(Sample-Reference)", 
                      row_names_gp = gpar(fontsize = row_fontsize),
                      top_annotation = ca,
                      cluster_columns = F,
                      cluster_rows = T)
# p_rtk

p_antiapop <- Heatmap(mat2plot[(rownames(mat2plot) %in% antiapop_genes) & (rowSums(!is.na(mat2plot)) >= num_nonna),],
                      row_title = "Pro-survival",
                      row_title_gp = gpar(fontsize = 12),
                      col = col_protein,
                      row_names_gp = gpar(fontsize = row_fontsize),
                      cluster_columns = F,
                      cluster_rows = T)
# p_antiapop

p_proapop <- Heatmap(mat2plot[(rownames(mat2plot) %in% proapop_genes) & (rowSums(!is.na(mat2plot)) >= num_nonna),],
                     row_title = "Pro-apoptosis",
                     row_title_gp = gpar(fontsize = 12),
                     col = col_protein,
                     row_names_gp = gpar(fontsize = row_fontsize),
                     cluster_columns = F,
                     cluster_rows = T)
# p_proapop

p_cdk <- Heatmap(mat2plot[(rownames(mat2plot) %in% cdk_genes) & (rowSums(!is.na(mat2plot)) >= num_nonna),],
                 row_title = "CDKs",
                 row_title_gp = gpar(fontsize = 12),
                 col = col_protein,
                 row_names_gp = gpar(fontsize = row_fontsize),
                 cluster_columns = F,
                 show_column_names = F,
                 cluster_rows = T)

p_mtor <- Heatmap(mat2plot[(rownames(mat2plot) %in% mtor_genes) & (rowSums(!is.na(mat2plot)) >= num_nonna),],
                 row_title = "mTORC1/2",
                 row_title_gp = gpar(fontsize = 12),
                 col = col_protein,
                 row_names_gp = gpar(fontsize = row_fontsize),
                 cluster_columns = F,
                 show_column_names = F,
                 cluster_rows = T)
# p_cdk
p_list <- p_rtk %v% p_mtor %v% p_antiapop %v% p_proapop %v% p_cdk 

## save heatmap to file
file2write <- paste0(dir_out, "Drug_Targets_Protein_Expression.", run_id, ".png")
png(file2write, width = 2000, height = 900, res = 150)
draw(p_list,
     row_title = "Druggable Target Family", row_title_gp = gpar(fontsize = 16),
     column_title = "Protein Abundance of the Druggable Targets", column_title_gp = gpar(fontsize = 16))
dev.off()


