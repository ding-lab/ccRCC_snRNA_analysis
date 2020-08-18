# Yige Wu @ WashU 2020 Feb
## plot a heatmap with CNV status for the frequently altered genes and their expression

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/load_data.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input RNA expression
rna_tab <- fread("./Resources/Bulk_Processed_Data/mRNA/RNA_rpkm_tumor_normal.tsv", data.table = F)
## input CNA matrix
cna_tab <- loadCNAstatus()
## input bulk meta data
bulk_meta_tab <- fread("./Resources/Bulk_Processed_Data/Meta_Data/cptac-metadata.csv")
## input metadata
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")
## input CNV genes
ccrcc_cna_genes_df <- readxl::read_excel(path = "./Resources/Knowledge/Known_Genetic_Alterations/Known_CNV.20200528.v1.xlsx", sheet = "Genes")
## set NonNA threshold
num_nonna <- 0

# set the aliquot IDs for bulk corresponding to the snRNA aliquots --------
case_ids2plot <- unique(idmetadata_df$Case[idmetadata_df$snRNA_available])
case_ids2plot

normal_bulk_aliquot_ids2plot <- mapvalues(x = case_ids2plot, from = bulk_meta_tab$Case.ID[bulk_meta_tab$Type == "Normal" & bulk_meta_tab$Set.A == "yes"], to = as.vector(bulk_meta_tab$RNA.ID[bulk_meta_tab$Type == "Normal" & bulk_meta_tab$Set.A == "yes"]))
normal_bulk_aliquot_ids2plot <- normal_bulk_aliquot_ids2plot[!(normal_bulk_aliquot_ids2plot %in% case_ids2plot)]
normal_bulk_aliquot_ids2plot

tumor_bulk_aliquot_ids2plot <- mapvalues(x = case_ids2plot, from = bulk_meta_tab$Case.ID[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"], to = as.vector(bulk_meta_tab$RNA.ID[bulk_meta_tab$Type == "Tumor" & bulk_meta_tab$Set.A == "yes"]))
tumor_bulk_aliquot_ids2plot

# plot gene expression data -------------------------------------
## specify the genes to be plotted 
genes2plot <- ccrcc_cna_genes_df$Gene_Symbol
genes2plot <- intersect(genes2plot, cna_tab$gene)
genes2plot
## filter the CNVs
cna_tab <- cna_tab[cna_tab$gene %in% genes2plot, c("gene", case_ids2plot)]
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


# make colors -------------------------------------------------------------
## get color corresponding to values
col_rna <- colorRamp2(c(quantile(mat2plot, 0.1, na.rm=T), quantile(mat2plot, 0.5, na.rm=T), quantile(mat2plot, 0.9, na.rm=T)),c("blue", "white", "red"))
# col_rna <- colorRamp2(c(-1.5, 0, 1.5),c("blue", "white", "red"))

# make column annotation --------------------------------------------------
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
bca = HeatmapAnnotation(df = col_anno_df, col = col_anno_colors, height = unit(1, "cm"))
tca = HeatmapAnnotation(Sample_Type = ifelse(grepl(x = colnames(mat2plot), pattern = "Tumor"), "Tumor", "Normal"),
                        col = list(Sample_Type = c("Tumor" = "red", "Normal" = "green")))

## plot heatmaps
p <- Heatmap(mat2plot[(rowSums(!is.na(mat2plot)) >= num_nonna),],
             width = unit(nrow(mat2plot), "cm"), height = unit(ncol(mat2plot)/2, "cm"),
             col = col_rna, na_col = "grey50",
             name = "log10(RPKM)", 
             row_names_gp = gpar(fontsize = 10),
             column_names_gp = gpar(fontsize = 10),
             bottom_annotation = bca,
             top_annotation = tca, 
             column_names_side = "top",
             show_heatmap_legend = F,
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
pdf(file2write, width = 20, height = 20)
print(p)
dev.off()
