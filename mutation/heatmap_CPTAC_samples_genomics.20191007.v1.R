# Yige Wu @ WashU 2019 Aug
## plot a heatmap with genomics data and proteomics data for given genes


###########################################
######## Source
###########################################
# source ------------------------------------------------------------------
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")

###########################################
######## FUNCTIONS
###########################################


###########################################
######## Set Variables
###########################################
# set parameters ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# set variables -----------------------------------------------------------
## plotting paramters
cap <- 3
breaks = seq(-(cap),cap, by=0.2)
## add color palette
color.palette <- colorRampPalette(rev(RColorBrewer::brewer.pal(10, "RdBu")))(length(breaks))

# plot MET, VEGFR  ---------------------------------------------------
mut_genes <- c(SMGs[["CCRCC"]])
cna_genes <- NULL
rna_genes <- c("VHL", "BAP1", "CA9",
               "HIF1A","EPAS1", 
               "VEGFA", "FGF2", "PDGFB", "FGFR1",
               "PDCD1", "CD274", "PDCD1LG2", "CTLA4")
pro_genes <- c("VHL", "BAP1", "CA9",
               "HIF1A","EPAS1", 
               "MET", "AXL",
               "VEGFA", ## VEGFA controlled by HIFs
               "FLT4", "FLT1", "KDR", "FLT3", 
               "FGF2", ## mostly FGF2 is reported
               "FGFR1",
               "PDGFA","PDGFB", ## PDGFB controlled by /HIFs
               "PDGFRA", "PDGFRB",
               "CDK4", "CDK6", "CDK7")
pho_genes_uniq <- c( "MET")
pho_genes <- pho_genes_uniq
rsds <- c("Y1234")
row_order <- c(paste0(rna_genes, "_RNA"),
               paste0(pro_genes, "_PRO"),
               paste0(pho_genes, "_", rsds))
row_order <- c("VHL_RNA", "VHL_PRO",
               "BAP1_RNA", "BAP1_PRO",
               "CA9_RNA", "CA9_PRO",
               "MET_PRO", "AXL_PRO",
               "HIF1A_RNA", "HIF1A_PRO",
               "EPAS1_RNA", "EPAS1_PRO",
               "VEGFA_RNA", "VEGFA_PRO",
               "FLT1_PRO",  "FLT4_PRO", "KDR_PRO",
               "FGF2_RNA", "FGF2_PRO",
               "FGFR1_RNA", "FGFR1_PRO",
               "PDGFB_RNA",
                "PDGFRB_PRO",
               "PDCD1_RNA", "CD274_RNA", "PDCD1LG2_RNA", "CTLA4_RNA",
               "CDK4_PRO", "CDK6_PRO", "CDK7_PRO")
fig_width <- 20
fig_height <- 8
nonNA_cutoff <- 0
version_tmp <- 1
if_cluster_row_tmp <- F
if_cluster_col_tmp <- F
if_col_name_tmp <- T

# plot Immune and metabolic genes ---------------------------------------------------
# mut_genes <- c(SMGs[["CCRCC"]])
# cna_genes <- c("")
# rna_genes <- c("PDGFRA", "PDCD1", "CD274", "PDCD1LG2", "CTLA4")
# pro_genes <- c("ENTPD1", "NT5E" , "POSTN","VEGFA", "PKM", "PRDX4", "FAP")
# pho_genes_uniq <- NULL
# pho_genes <- pho_genes_uniq
# rsds <- NULL
# row_order <- c(paste0(rna_genes, "_RNA"), paste0(pro_genes, "_PRO"))
# fig_width <- 20
# fig_height <- 5
# nonNA_cutoff <- 20
# version_tmp <- 1
# if_cluster_row_tmp <- T
# if_cluster_col_tmp <- T
# if_col_name_tmp <- T

###########################################
######## Plotting
###########################################

# bussiness ------------------------------------------------------------------
geneA <- paste(head(unique(c(mut_genes, cna_genes)), 5), collapse = "_")
geneB <- paste(head(unique(c(rna_genes, pro_genes, pho_genes)), 8), collapse = "_")
phosphosite <- paste0(rsds, collapse = "_")

ann_colors <- list()

# input data first because different for each cancer type --------------------------------------------------------------
## input mutation matrix
maf <- loadMaf()
somatic_mat <- generate_somatic_mutation_matrix(pair_tab = mut_genes, maf = maf)

## mutation needs to show both geneA and geneB
somatic_mat <- somatic_mat[somatic_mat$Hugo_Symbol %in% mut_genes,]

## input CNA matrix
cna_tab <- loadCNAstatus()
cna_tab <- cna_tab[cna_tab$gene %in% cna_genes, ]


## load RNA
rna_tab <- loadRNA()
rna_tab <- rna_tab[rna_tab$gene %in% rna_genes,]

## input protein data
pro_tab <- loadParseProteomicsData(expression_type  = "PRO", sample_type = "tumor")
pho_tab <- loadParseProteomicsData(expression_type  = "PHO", sample_type = "tumor")

pro_tab <- pro_tab[pro_tab$Gene %in% pro_genes,]
pho_tab <- pho_tab[pho_tab$Gene %in% pho_genes & pho_tab$Phosphosite %in% rsds,]

# make the annotation columns for each sample -----------------------------
partIDs <- colnames(pho_tab)[!(colnames(pho_tab) %in% c("Gene", "Phosphosite", "Peptide_ID"))]
col_anno <- data.frame(partID = partIDs)

if (nrow(somatic_mat) > 0){
  somatic_mat.m <- melt(somatic_mat, id.vars = "Hugo_Symbol")
  somatic_mat.m %>% head()
  colnames(somatic_mat.m) <- c("Gene", "partID", "variant_class")
  
  ## distinguish by missense and truncation
  somatic_mat.m$variant_class[is.na(somatic_mat.m$variant_class)] <- ""
  somatic_mat.m$variant_class_sim <- "other_mutation"
  somatic_mat.m$variant_class_sim[somatic_mat.m$variant_class == ""] <- "wild_type"
  somatic_mat.m$variant_class_sim[somatic_mat.m$variant_class  == "Silent"] <- "silent"
  somatic_mat.m$variant_class_sim[grepl(x = somatic_mat.m$variant_class, pattern = "Missense_Mutation")] <- "missense"
  somatic_mat.m$variant_class_sim[grepl(x = somatic_mat.m$variant_class, pattern = "Nonsense_Mutation|Frame_Shift_Ins|Frame_Shift_Del")] <- "truncation"
  somatic_mat.m$variant_class_sim[sapply(X = somatic_mat.m$variant_class, FUN = function(v) (grepl(pattern = "Nonsense_Mutation|Frame_Shift_Ins|Frame_Shift_Del", x = v) & grepl(pattern = "Missense_Mutation", x = v)))] <- "missense&truncation"
  
  for (gene in unique(somatic_mat.m$Gene[somatic_mat.m$variant_class_sim != "wild_type"])) {
    somatic_mat2merge <- somatic_mat.m[somatic_mat.m$Gene == gene, c("partID", "variant_class_sim")]
    colnames(somatic_mat2merge) <- c("partID", paste0("somatic.", gene))
    col_anno <- merge(col_anno, somatic_mat2merge, by = c("partID"), all.x = T)
  }
} else {
  print("no mutation!")
}

## CNA needs to show both geneA and geneB; levels: amplification, deletion, neutral
if (nrow(cna_tab) > 0) {
  cna_tab.m <- melt(cna_tab, id.vars = "gene")
  colnames(cna_tab.m) <- c("gene", "partID", "CNA")
  
  for (gene in intersect(cna_genes, unique(cna_tab.m$gene[cna_tab.m$CNA != "neutral"]))) {
    cna_mat2merge <- cna_tab.m[cna_tab.m$gene == gene, c("partID", "CNA")]
    colnames(cna_mat2merge) <- c("partID", paste0("CNA.", gene))
    col_anno <- merge(col_anno, cna_mat2merge, by = c("partID"), all.x = T)
  }
} else {
  print("no CNA!")
}

# input sample mapping file -----------------------------------------------
cptac_sample_map <- fread(input = "./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC_expression_matrices/cptac-metadata.csv", data.table = F)
cptac_sample_map <- data.frame(cptac_sample_map)
specimen2case_map <- cptac_sample_map %>%
  dplyr::filter(Type == "Tumor") %>%
  select(Case.ID, Specimen.Label)
rownames(specimen2case_map) <- specimen2case_map$Specimen.Label


# input xcell result and add immnue group ------------------------------------------------------
xcell_tab <- readxl::read_excel("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC Manuscript/CPTAC3-ccRCC-SupplementaryTables_Final/Table S7.xlsx", sheet = "xCell Signatures", skip = 2)
## immnue group 
immune_groups <- unlist(xcell_tab[1,2:ncol(xcell_tab)])
immune_group_df <- data.frame(immune_group = immune_groups, partID = specimen2case_map[names(immune_groups), "Case.ID"])
col_anno <- merge(col_anno, immune_group_df, by = c("partID"), all.x = T)

# add info for original tumor sample at WashU -------------------------------------------------------
shipped_sample_info <- read_excel("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/Sample_Info/02_Post_request_CPTAC3_shipment/JHU Distribution Manifest_ccRCC material to WUSTL_6.17.2019_original_segment.xlsx")
col_anno$Tumor_Original_at_WashU <- as.character(col_anno$partID %in% shipped_sample_info$`Subject ID`)


col_anno %>% head()

for (gene in cna_genes) {
  if (paste0("CNA.", gene) %in% colnames(col_anno)) {
    col_anno <- col_anno[order(col_anno[, paste0("CNA.", gene)], decreasing = T),]
    ann_colors[[paste0("CNA.", gene)]] <-  c(amplification = "#E41A1C", deletion = "#377EB8", "neutral" = "grey")
  }
}
for (gene in mut_genes) {
  if (paste0("somatic.", gene) %in% colnames(col_anno)) {
    col_anno <- col_anno[order(col_anno[, paste0("somatic.", gene)], decreasing = T),]
    ann_colors[[paste0("somatic.", gene)]] <- c(missense = "#E41A1C", truncation = "#377EB8", wild_type = "white", "missense&truncation" = "#6A3D9A", other_mutation = "#FF7F00", silent = "#33A02C")
  }
}

# if ("MET_PRO" %in% colnames(col_anno)) {
#   ann_colors[["MET_Y1234_PHO"]] <- c("TRUE" = "#E41A1C", "FALSE" = "grey", "NA" = "grey")
#   ann_colors[["MET_PRO"]] <- c("TRUE" = "#E41A1C", "FALSE" = "grey", "NA" = "grey")
# }

if ("immune_group" %in% colnames(col_anno)) {
  ann_colors[["immune_group"]] <- c("CD8+ inflamed" = "#E41A1C", "CD8- inflamed" = "#FF7F00", "VEGF immune-desert" = "#377EB8", "Metabolic immune-desert" = "#6A3D9A", "NA" = "white")
}

if ("Tumor_Original_at_WashU" %in% colnames(col_anno)) {
  ann_colors[["Tumor_Original_at_WashU"]] <-  c("TRUE" = "#E41A1C", "FALSE" = "grey", "NA" = "grey")
}

# col_anno <- col_anno[order(col_anno$variant_class_sim),]
col_anno %>% head()
rownames(col_anno) <- col_anno$partID

# make the matrix of values showing in heatmap ----------------------------
sup_tab_can <- NULL

if (nrow(rna_tab) > 0) {
  rna_tab.m <- melt(rna_tab, id.vars = "gene")
  colnames(rna_tab.m) <- c("Gene", "partID", "exp_value")
  rna_tab.m$Phosphosite <- "RNA"
  sup_tab_can <- rbind(sup_tab_can, rna_tab.m[,c("Gene", "Phosphosite", "partID", "exp_value")])
}

if (nrow(pro_tab) > 0) {
  pro_tab.m <- melt(pro_tab, id.vars = "Gene")
  pro_tab.m %>% head()
  colnames(pro_tab.m) <- c("Gene", "partID", "exp_value")
  pro_tab.m$Phosphosite <- "PRO"
  sup_tab_can <- rbind(sup_tab_can, pro_tab.m[,c("Gene", "Phosphosite", "partID", "exp_value")])
}

if (nrow(pho_tab) > 0) {
  pho_tab <- pho_tab[!duplicated(pho_tab[, c("Gene", "Phosphosite")]),!(colnames(pho_tab) %in% c("Peptide_ID"))]
  pho_tab.m <- melt(pho_tab, id.vars = c("Gene", "Phosphosite"))
  pho_tab.m %>% head()
  colnames(pho_tab.m) <- c("Gene", "Phosphosite", "partID", "exp_value")
  sup_tab_can <- rbind(sup_tab_can, pho_tab.m[,c("Gene", "Phosphosite", "partID", "exp_value")])
}

sup_tab_can$id_row <- paste0(sup_tab_can$Gene, "_", sup_tab_can$Phosphosite)
sup_tab_can$exp_value <- as.numeric(as.vector(sup_tab_can$exp_value))
sup_tab_can <- unique(sup_tab_can)

## make the matrix for the heatmap body
df_value <- dcast(data = sup_tab_can, id_row ~ partID, value.var = "exp_value")

df_value %>% head()
mat_value <- as.matrix(df_value[,-1])
rownames(mat_value) <- df_value$id_row
head(mat_value)
# mat_value <- mat_value[,colSums(!is.na(mat_value)) > 0]
mat_value %>% head()

## order the matrix column
partIDs_ordered <- intersect(as.vector(rownames(col_anno)), colnames(mat_value))
col_anno$partID <- NULL
mat_value <- mat_value[, partIDs_ordered]

## order the matrix rows
# row_order <- c(paste0(rep(unique(c(rna_genes, pro_genes, phog_genes)), 3), "_", rep(c("RNA", "PRO", "collapsed_PHO"), length(unique(c(rna_genes, pro_genes, pho_genes))) + c(0,0,0))), paste0(pho_genes, "_", rsds))
if (length(row_order) > 1) {
  mat_value <- mat_value[intersect(row_order, rownames(mat_value)),]
  mat_value <- mat_value[rowSums(!is.na(mat_value)) >= nonNA_cutoff, ]
  mat_value <- mat_value[,colSums(!is.na(mat_value)) >= 1]
} else {
  mat_value <- matrix(data = mat_value, nrow = 1, dimnames = list(row_order, names(mat_value)))
}

# plotting ----------------------------------------------------------------
my_heatmap <- pheatmap(mat_value, 
                       color = color.palette,
                       annotation_col = col_anno,
                       scale = "row",
                       na_col = "white",
                       show_colnames = if_col_name_tmp,
                       cluster_rows=if_cluster_row_tmp, 
                       cluster_cols=if_cluster_col_tmp, 
                       gaps_row = c(6, 8, 12, 17, 21, 23, 28),
                       annotation_colors = ann_colors)

fn <- paste0(dir_out, "CPTAC_ccRCC_Genomics_Landscape", ".", run_id, ".png")
save_pheatmap_png(x = my_heatmap, 
                  filename = fn, 
                  width = 2500, height = 1200)


fn <- paste0(dir_out, "CPTAC_ccRCC_Genomics_Landscape", ".", run_id, ".pdf")
save_pheatmap_pdf(x = my_heatmap, 
                  filename = fn, 
                  width = fig_width, height = fig_height)




