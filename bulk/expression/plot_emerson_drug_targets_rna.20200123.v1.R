# Yige Wu @ WashU 2019 Nov
## plot a heatmap with proteomics data from the discovery set data freeze

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")


# set run id --------------------------------------------------------------
version_tmp <- 2
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)


# set output directory ----------------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input protein data ------------------------------------------------------
rna_tab <- fread("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/mRNA/RNA_rpkm_tumor_normal.tsv", data.table = F)

# input bulk meta data ----------------------------------------------------
bulk_meta_tab <- fread("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/cptac-metadata.csv")

# get the aliquot IDs for bulk corresponding to the snRNA aliquots --------
normal_bulk_aliquot_ids2plot <- bulk_meta_tab$RNA.ID[bulk_meta_tab$Set.A == "yes" & bulk_meta_tab$Type == "Normal" & !is.na(bulk_meta_tab$RNA.ID)]
normal_bulk_aliquot_ids2plot

case_ids2plot <- mapvalues(x = normal_bulk_aliquot_ids2plot, from = bulk_meta_tab$RNA.ID, to = bulk_meta_tab$Case.ID)
case_ids2plot

tumor_bulk_aliquot_ids2plot <- mapvalues(x = case_ids2plot, from = bulk_meta_tab$Case.ID[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"], to = as.vector(bulk_meta_tab$RNA.ID[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"]))
tumor_bulk_aliquot_ids2plot


# set NonNA threshold -----------------------------------------------------
num_nonna <- 40
row_fontsize <- 9

# input gene list ---------------------------------------------------------
## targets by the IAP inhibitors: https://www.medchemexpress.com/Targets/IAP.html?src=google-product&gclid=EAIaIQobChMIr9eU1IWY5wIVEPDACh1IGA-_EAAYAiAAEgKxLPD_BwE
# IAP (Inhibitors of Apoptosis) is a family of functionally and structurally related proteins, which serve as endogenous inhibitors of programmed cell death (apoptosis). A common feature of all IAPs is the presence of a BIR in one to three copies. The human IAP family consists of 8 members, and IAP homologs have been identified in numerous organisms. The members of the IAPs included IAPs, Cp-IAP, Op-IAP, XIAP, c-IAPl, C-IAP2, NAIP, Livin and Survivin. The best characterized IAP is XIAP, which binds caspase-9, caspase-3 and caspase 7, thereby inhibiting their activation and preventing apoptosis. Also cIAP1 and cIAP2 have been shown to bind caspases, although how the IAPs inhibit apoptosis mechanistically at the molecular level is not completely understood.
IAP_genes <- c("BIRC2", "BIRC3", "XIAP", "BIRC7", "BIRC5")
cdk_genes <- c("CDK1", "CDK2", "CDK4", "CDK5", "CDK6", "CDK7", "CDK9", "WEE1")
mtor_genes <- c("MTOR", "AKTS1", "DEPTOR", "RPTOR", "MLST8", "MAPKAP1", "PRR5", "RICTOR")
## targets by carbozantinib: https://www.drugbank.ca/drugs/DB08875
rtk_genes <- c("MET", "AXL", "FLT1", "KDR", "FLT3", "KIT", "RET", "NTRK2", "TEK")

# plot gene expression data -------------------------------------
## specify the genes to be plotted 
genes2plot <- c(IAP_genes, cdk_genes, mtor_genes, rtk_genes)

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
## make column annotation
ca = HeatmapAnnotation(Sample_Type = ifelse(colnames(mat2plot) %in% tumor_bulk_aliquot_ids2plot, "Tumor", "Normal"),
                       col = list(Sample_Type = c("Tumor" = "red", "Normal" = "green")))

## get color corresponding to values
col_rna <- colorRamp2(c(min(mat2plot, na.rm=T), mean(mat2plot, na.rm=T), max(mat2plot, na.rm=T)),c("blue", "white", "red"))
col_rna <- colorRamp2(c(quantile(mat2plot, 0.1, na.rm=T), mean(mat2plot, na.rm=T), quantile(mat2plot, 0.9, na.rm=T)),c("blue", "white", "red"))
col_rna <- colorRamp2(c(-1.5, 0, 1.5),c("blue", "white", "red"))
col_rna <- colorRamp2(c(-2, 0, 2),c("blue", "white", "red"))

## plot heatmaps
p_rtk <- Heatmap(mat2plot[(rownames(mat2plot) %in% rtk_genes) & (rowSums(!is.na(mat2plot)) >= num_nonna),],
                      row_title = "RTKs",
                      row_title_gp = gpar(fontsize = 12),
                      col = col_rna,
                      name = "log10(RPKM)", 
                      row_names_gp = gpar(fontsize = row_fontsize),
                      top_annotation = ca,
                      cluster_columns = F,
                      cluster_rows = T)

p_IAP <- Heatmap(mat2plot[(rownames(mat2plot) %in% IAP_genes) & (rowSums(!is.na(mat2plot)) >= num_nonna),],
                      row_title = "IAPs",
                      row_title_gp = gpar(fontsize = 12),
                      col = col_rna,
                      row_names_gp = gpar(fontsize = row_fontsize),
                      cluster_columns = F,
                      cluster_rows = T)

p_cdk <- Heatmap(mat2plot[(rownames(mat2plot) %in% cdk_genes) & (rowSums(!is.na(mat2plot)) >= num_nonna),],
                 row_title = "CDKs",
                 row_title_gp = gpar(fontsize = 12),
                 col = col_rna,
                 row_names_gp = gpar(fontsize = row_fontsize),
                 cluster_columns = F,
                 show_column_names = F,
                 cluster_rows = T)

p_mtor <- Heatmap(mat2plot[(rownames(mat2plot) %in% mtor_genes) & (rowSums(!is.na(mat2plot)) >= num_nonna),],
                 row_title = "mTORC1/2",
                 row_title_gp = gpar(fontsize = 12),
                 col = col_rna,
                 row_names_gp = gpar(fontsize = row_fontsize),
                 cluster_columns = F,
                 show_column_names = F,
                 cluster_rows = T)
# p_cdk
p_list <- p_rtk %v% p_mtor %v% p_IAP %v% p_cdk 

## save heatmap to PDF
file2write <- paste0(dir_out, "Drug_Targets_mRNA_Expression_log10RPKM.", run_id, ".pdf")
pdf(file2write, width = 20, height = 5)
draw(p_list,
     row_title = "Druggable Target Family", row_title_gp = gpar(fontsize = 16),
     column_title = "mRNA Expression of the Druggable Targets", column_title_gp = gpar(fontsize = 16))
dev.off()


## save heatmap to PNG
file2write <- paste0(dir_out, "Drug_Targets_mRNA_Expression_log10RPKM.", run_id, ".png")
png(file2write, width = 2000, height = 600, res = 150)
draw(p_list,
     row_title = "Druggable Target Family", row_title_gp = gpar(fontsize = 16),
     column_title = "mRNA Expression of the Druggable Targets", column_title_gp = gpar(fontsize = 16))
dev.off()


