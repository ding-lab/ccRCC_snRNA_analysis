# Yige Wu @WashU Oct 2019
## for plotting the scaled snRNA expression with bulk protein to show correlation across cell clusters

# source ------------------------------------------------------------------
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")

# set run id  ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)


# set output directory ----------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# input integrated data ---------------------------------------------------
object2plot <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/integrate_seurat_objects/20191021.v1/Renal_Integrated.20191021.v1.RDS")
DefaultAssay(object2plot) <- "RNA"
aliquot_ids <- unique(object2plot@meta.data$orig.ident)


# get number of barcodes per cluster --------------------------------------
sn_cell_num_tab <- data.frame(object2plot@meta.data %>%
                                select(orig.ident, seurat_clusters) %>%
                                table())
sn_cell_num_tab %>%
  head()
colnames(sn_cell_num_tab) <- c("snRNA_Aliquot_ID", "Cluster", "Num_Cluster_Barcode")
sn_cell_sum_tab <- sn_cell_num_tab %>%
  group_by(snRNA_Aliquot_ID) %>%
  summarise(Num_Aliquot_Barcode = sum(Num_Cluster_Barcode))
sn_cell_num_tab <- merge(sn_cell_num_tab, sn_cell_sum_tab, by = c("snRNA_Aliquot_ID"), all.x = T)
sn_cell_num_tab <- sn_cell_num_tab %>%
  mutate(Perc_Cluster_Barcode = Num_Cluster_Barcode/Num_Aliquot_Barcode)

# input cluster cell type assignment --------------------------------------
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - AllCluster2Cell_Type.20191022.v1.tsv", data.table = F)
nonimmune_clusters <- cluster2celltype_tab %>%
  filter(Is_Immune == "No") %>%
  select(Cluster)
nonimmune_clusters <- nonimmune_clusters$Cluster

# input bulk protein data ------------------------------------------------------
pro_tab <- loadParseProteomicsData(expression_type  = "PRO", sample_type = "tumor")
pro_mat <- pro_tab
rownames(pro_mat) <- pro_tab$Gene
pro_mat$Gene <- NULL
pro_mat <- as.matrix(pro_mat)

# input bulk RNA data -----------------------------------------------------
rna_tab <- loadRNA()
rna_mat <- rna_tab
rownames(rna_mat) <- rna_tab$gene
rna_mat$gene <- NULL
rna_mat <- as.matrix(rna_mat)

# input meta data ---------------------------------------------------------
meta_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/make_meta_data/meta_data.20190924.v1.tsv", data.table = F)
meta_tab_in_use <- meta_tab %>%
  dplyr::filter(meta_tab$Specimen.ID.snRNA %in% aliquot_ids)
rownames(meta_tab_in_use) <- meta_tab_in_use$Case.ID

# input xcell result and add immnue group ------------------------------------------------------
xcell_tab <- readxl::read_excel("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC Manuscript/CPTAC3-ccRCC-SupplementaryTables_Final/Table S7.xlsx", sheet = "xCell Signatures", skip = 2)

## immnue group 
immune_groups <- unlist(xcell_tab[1,2:ncol(xcell_tab)])
immune_groups <- immune_groups[meta_tab_in_use$Specimen.ID.bulk]
immune_groups
names(immune_groups) <- meta_tab_in_use$Specimen.ID.snRNA

# input selected clinical data --------------------------------------------
clinical_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/generate_clinical_table/20191017.v1/snRNA_ccRCC_Clinicl_Table.20191017.v1.tsv", data.table = F)
clinical_tab <- clinical_tab %>%
  filter(Case %in% meta_tab_in_use$Case.ID)
clinical_tab <- merge(meta_tab_in_use, clinical_tab, by.x = c("Case.ID"), by.y = c("Case"), all.x = T)
clinical_tab_outcome_initial_mat <- str_split_fixed(string = clinical_tab$Outcome_of_Initial_Treatment, pattern = "\\|", n = 3)
clinical_tab_outcome_followup_mat <- str_split_fixed(string = clinical_tab$Outcome_at_Followup_Completion, pattern = "\\|", n = 3)
clinical_tab_outcome_followup_mat[which(clinical_tab_outcome_followup_mat == "Not Applicable", arr.ind = T)] <- clinical_tab_outcome_initial_mat[which(clinical_tab_outcome_followup_mat == "Not Applicable", arr.ind = T)]
clinical_tab_outcome_followup_mat
colnames(clinical_tab_outcome_followup_mat) <- paste0("Outcome_", seq(from = 12, to = 36, by = 12), "Months")
clinical_tab <- cbind(clinical_tab, clinical_tab_outcome_followup_mat)

# reformat bulk protein and RNA data colnames to aliquot ids ---------------------------
pro_mat <- pro_mat[,meta_tab_in_use$Case.ID]
pro_mat %>% head()
colnames(pro_mat) <- meta_tab_in_use$Specimen.ID.snRNA

rna_mat <- rna_mat[,meta_tab_in_use$Case.ID]
rna_mat %>% head()
colnames(rna_mat) <- meta_tab_in_use$Specimen.ID.snRNA

# set genes to plot -------------------------------------------------------
gene_tmp <- "CD274"
genes2plot <- c("PDCD1", "CD274", "PDCD1LG2", "CTLA4", "CD38",
                "OXPHOS", "PRDX4", "PKM",
                "PDGFRA", "FAP", "ENTPD1", "NT5E", "POSTN",
                "CTNNB1", "RAP1")
genes2plot <- c("MET", "CA9")
genes2plot <- SMGs[["CCRCC"]]

for (gene_tmp in genes2plot) {
  # gather snRNA data as well as protein data into a dataframe -------------------------------
  tab2p <- NULL
  for (aliquot_id_tmp in names(immune_groups)) {
    object2plot_tmp <- subset(x = object2plot,subset = orig.ident == aliquot_id_tmp)
    p <- DotPlot(object = object2plot_tmp, features = gene_tmp, group.by = c("seurat_clusters"))
    p$data$cluster_cell_type <- plyr::mapvalues(p$data$id, from = cluster2celltype_tab$Cluster, to = cluster2celltype_tab$Enriched_Cell_Type_Abbr)
    p$data$is_immune <- ifelse(p$data$id %in% nonimmune_clusters, "Non_Immune", "Immune")
    p$data$is_immune <- factor(p$data$is_immune, levels = c("Non_Immune", "Immune"))
    p$data$aliquot_id <- aliquot_id_tmp
    tab2p <- rbind(tab2p, p$data)
  }
  tab2p <- merge(tab2p, sn_cell_num_tab[, c("snRNA_Aliquot_ID", "Perc_Cluster_Barcode", "Cluster")], 
                 by.x = c("aliquot_id", "id"), by.y = c("snRNA_Aliquot_ID", "Cluster"), all.x = T)
  tab2p <- tab2p %>%
    mutate(pct.exp.time.pct.cluster = (pct.exp/100)*(Perc_Cluster_Barcode))
  
  # make the matrix with values for color -----------------------------------
  mat_color <- dcast(data = tab2p, formula = aliquot_id ~ cluster_cell_type, value.var = "avg.exp")
  rownames(mat_color) <- mat_color$aliquot_id
  mat_color$aliquot_id <- NULL
  mat_color <- as.matrix(mat_color)
  
  # make the matrix with values for circle size -----------------------------------
  mat_cluster_size <- dcast(data = tab2p, formula = aliquot_id ~ cluster_cell_type, value.var = "Perc_Cluster_Barcode")
  rownames(mat_cluster_size) <- mat_cluster_size$aliquot_id
  mat_cluster_size$aliquot_id <- NULL
  mat_cluster_size <- as.matrix(mat_cluster_size)
  
  mat_size <- dcast(data = tab2p, formula = aliquot_id ~ cluster_cell_type, value.var = "pct.exp.time.pct.cluster")
  rownames(mat_size) <- mat_size$aliquot_id
  mat_size$aliquot_id <- NULL
  mat_size <- as.matrix(mat_size)
  
  # make color transformation function --------------------------------------
  summary(tab2p$avg.exp)
  col_fun = circlize::colorRamp2(seq(from = 0, to = max(tab2p$avg.exp), length.out = 4), RColorBrewer::brewer.pal(n = 4, name = "YlOrRd"))
  
  # make row annotation ------------------------------------------------
  ## meke the color transformation function for bulk protein
  protein_col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  
  ## meke the color transformation function for bulk mRNA
  bulk_mrna_col_fun = colorRamp2(0:8, RColorBrewer::brewer.pal(n = 9, name = "YlOrRd"))
  
  ## make the row annotation
  if (gene_tmp %in% rownames(pro_mat)) {
    ra = rowAnnotation(Bulk_Protein = pro_mat[gene_tmp, rownames(mat_color)], 
                       Bulk_mRNA = rna_mat[gene_tmp, rownames(mat_color)],
                       Immune_Subtype = immune_groups,
                       Outcome_12Mon = mapvalues(x = rownames(mat_color), from = clinical_tab$Specimen.ID.snRNA, to = as.vector(clinical_tab$Outcome_12Months)),
                       Outcome_24Mon = mapvalues(x = rownames(mat_color), from = clinical_tab$Specimen.ID.snRNA, to = as.vector(clinical_tab$Outcome_24Months)),
                       Outcome_36Mon = mapvalues(x = rownames(mat_color), from = clinical_tab$Specimen.ID.snRNA, to = as.vector(clinical_tab$Outcome_36Months)),
                       col = list(Bulk_Protein = protein_col_fun,
                                  Bulk_mRNA = bulk_mrna_col_fun,
                                  Immune_Subtype = c("CD8+ inflamed" = "#E41A1C", "CD8- inflamed" = "#377EB8",
                                                     "Metabolic immune-desert" = "#33A02C", "VEGF immune-desert" = "#FDBF6F"),
                                  Outcome_12Mon = c("Complete Remission" = "#7fc97f", "Patient Deceased" = "black", "Unknown" = "grey50", "Not Applicable" = "white"),
                                  Outcome_24Mon = c("Complete Remission" = "#7fc97f", "Patient Deceased" = "black", "Unknown" = "grey50", "Not Applicable" = "white"),
                                  Outcome_36Mon = c("Complete Remission" = "#7fc97f", "Patient Deceased" = "black", "Unknown" = "grey50", "Not Applicable" = "white")))
  } else {
    ra = rowAnnotation(Bulk_mRNA = rna_mat[gene_tmp, rownames(mat_color)],
                       Immune_Subtype = immune_groups,
                       col = list(Bulk_mRNA = bulk_mrna_col_fun,
                                  Immune_Subtype = c("CD8+ inflamed" = "#E41A1C", "CD8- inflamed" = "#377EB8",
                                                     "Metabolic immune-desert" = "#33A02C", "VEGF immune-desert" = "#FDBF6F")))
  }

  
  # make column annotation ------------------------------------------------
  ca = HeatmapAnnotation(Is_Immune_Cluster = (colnames(mat_color) %in% cluster2celltype_tab$Enriched_Cell_Type_Abbr[cluster2celltype_tab$Is_Immune == "Yes"]),
                         col = list(Is_Immune_Cluster = c("TRUE" = "black", "FALSE" = "white")))
  
  
  # set column order --------------------------------------------------------
  column_order_vec <- cluster2celltype_tab$Is_Immune
  names(column_order_vec) <- cluster2celltype_tab$Enriched_Cell_Type_Abbr
  column_order_vec <- column_order_vec[colnames(mat_color)]
  
  
  # set row order -----------------------------------------------------------
  row_order_vec <- immune_groups
  
  # plot dot plot per sample per gene after cell type assignment-----------------------------------------------------------
  p <- Heatmap(mat_color, 
               name = "avg.exp", 
               col = col_fun, 
               rect_gp = gpar(type = "none"), 
               cell_fun = function(j, i, x, y, width, height, fill) {
                 # grid.rect(x = x, y = y, width = width, height = height,
                 #           gp = gpar(col = "grey", fill = NA))
                 grid.circle(x = x, y = y, r = sqrt(mat_cluster_size[i, j])*2.5* min(unit.c(width, height)), 
                             gp = gpar(col = "black", fill = NA))
                 grid.circle(x = x, y = y, r = sqrt(mat_size[i, j])*2.5 * min(unit.c(width, height)), 
                             gp = gpar(fill = col_fun(mat_color[i, j]), col = NA))
                 # if (mat_size[i, j] > 1) {
                 #   grid.text(sprintf("%.1f", 100*mat_size[i, j]), x, y, gp = gpar(fontsize = 10))
                 # }
               }, 
               row_names_side = "left", 
               left_annotation = ra,
               row_order = order(immune_groups),
               column_title = paste0(gene_tmp, " Single Nuclei Expression Per Cell Cluster"),
               top_annotation = ca,
               column_order = order(column_order_vec),
               cluster_rows = FALSE, 
               cluster_columns = FALSE,
               show_row_names = T, 
               show_column_names = T)
  
  
  file2write <- paste0(dir_out, gene_tmp, "_Exp_Per_Cluster_w_Bulk_Protein_mRNA.", run_id, ".png")
  png(file2write, width = 3000, height = 1000, res = 150)
  print(p)
  dev.off()
}
