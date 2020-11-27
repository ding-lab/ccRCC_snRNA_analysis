# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the average expression calculated (SCT)
avgexp_df <- fread(input = "./Resources/Analysis_Results/average_expression/averageexpression_sct_usescale_byindividualcluster_bycellgroup7_byaliquot_on_katmai/20200917.v1/avgexp.SCT.bycellgroup.byaliquot.bycluster.31_aliquot_integration.20200917.v1.tsv", data.table = F)
## input the barcode-cell-type table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups/20200917.v2/31Aliquot.Barcode2CellType.20200917.v2.tsv", data.table = F)
## barcode 2 individual cluster id
barcode2cluster_df <- fread(data.table = F, input = "./Resources/Analysis_Results/fetch_data/fetch_data_by_individual_sample/20200717.v1/Barcode2MetaData.20200717.v1.tsv")
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")

# specify genes to filter -------------------------------------------------
## input kidney-specific EMT genes
emt_genes_df <- fread(data.table = F, input = "./Resources/Analysis_Results/dependencies/combine_pt_with_emt_markers_all/20200920.v1/Kidney_Specific_EMT_Genes.20200920.v1.tsv")
## add name for the marker groups
emt_genes_df <- emt_genes_df %>%
  mutate(Text_Gene_Group = ifelse(Gene_Group2 == "Tumor cells", 
                                  "Tumor-cell\nmarkers", 
                                  ifelse(Gene_Group2 %in% c("Epithelial", "Proximal tubule"),
                                         "Epithelial/\nproximal-tubule\nmarkers", paste0(Gene_Group2, "\nmarkers"))))

genes2filter <- emt_genes_df$Gene

# count cell number and filter clusters -----------------------------------
barcode2celltype_df <- merge(barcode2celltype_df, barcode2cluster_df, by.x = c("orig.ident", "individual_barcode"), by.y = c("aliquot", "individual_barcode"), all.x = T)
barcode2celltype_df <- barcode2celltype_df %>%
  mutate(id_bycluster_bycellgroup_byaliquot = paste0(orig.ident, "_", seurat_cluster_id, "_",Cell_group7))
cellcount_bycluster_df <- barcode2celltype_df %>%
  select(id_bycluster_bycellgroup_byaliquot) %>%
  table() %>%
  as.data.frame() %>%
  rename(id_bycluster_bycellgroup_byaliquot_original = ".") %>%
  mutate(id_bycluster_bycellgroup_byaliquot = gsub(x = id_bycluster_bycellgroup_byaliquot_original, pattern = " |\\-", replacement = "."))

# format expression data --------------------------------------------------
plot_data_long_df <- avgexp_df %>%
  filter(V1 %in% genes2filter) %>%
  melt() %>%
  mutate(id_bycluster_bycellgroup_byaliquot = gsub(x = variable, pattern = "SCT.", replacement = "")) %>%
  mutate(aliquot = str_split_fixed(string = id_bycluster_bycellgroup_byaliquot, pattern = "_", n = 3)[,1]) %>%
  mutate(id_cluster = str_split_fixed(string = id_bycluster_bycellgroup_byaliquot, pattern = "_", n = 3)[,2]) %>%
  mutate(cellgroup = str_split_fixed(string = id_bycluster_bycellgroup_byaliquot, pattern = "_", n = 3)[,3])
plot_data_long_df$Cell_count <- mapvalues(x = plot_data_long_df$id_bycluster_bycellgroup_byaliquot, from = cellcount_bycluster_df$id_bycluster_bycellgroup_byaliquot, to = as.vector(cellcount_bycluster_df$Freq))
plot_data_long_df$Cell_count <- as.numeric(as.vector(plot_data_long_df$Cell_count))
plot_data_long_df <- plot_data_long_df %>%
  dplyr::filter(Cell_count >= 30) %>%
  dplyr::filter(cellgroup %in% c("Tumor.cells", "Transitional.cells", "Tumor.like.cells"))
plot_data_long_df$id_aliquot_wu <- mapvalues(x = plot_data_long_df$aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
plot_data_long_df <- plot_data_long_df %>%
  dplyr::mutate(id_bycluster_bycellgroup_byaliquot_new = paste0(id_aliquot_wu, "_", id_cluster, "_", cellgroup))
## make matrix
plot_data_wide_df <- dcast(data = plot_data_long_df, formula = V1 ~ id_bycluster_bycellgroup_byaliquot_new, value.var = "value")
plot_data_mat <- as.matrix(plot_data_wide_df[,-1])
rownames(plot_data_mat) <- plot_data_wide_df$V1

# filter genes based on variation -----------------------------------------
sd_bygene_df <- data.frame(SD = apply(plot_data_mat,1, sd, na.rm = TRUE), gene = rownames(plot_data_mat))
sd_bygene_df$Cell_Group2 <- mapvalues(x = sd_bygene_df$gene, from = emt_genes_df$Gene, to = as.vector(emt_genes_df$Gene_Group2))
sd_bygene_df <- sd_bygene_df %>%
  arrange(desc(SD))
genes_plot_mesenchymal <- head(as.vector(sd_bygene_df$gene[sd_bygene_df$Cell_Group2 == "Mesenchymal"]), n = 5)
genes_plot_mesenchymal
sd_bygene_epithelial_df <- sd_bygene_df %>%
  filter(Cell_Group2 %in% c("Proximal tubule")) %>%
  arrange(desc(SD))
genes_plot_epithelial <- head(x = as.vector(sd_bygene_epithelial_df$gene), n = 5)
genes_plot_other <- genes_plot_epithelial
genes_plot <- c(genes_plot_mesenchymal, genes_plot_other)
plot_data_mat <- plot_data_mat[genes_plot,]
## make mesenchymal score
genes_mesenchymal_score <- c("FN1", "CDH2", "VIM", "FOXC2", "SNAI2")
scores_mesenchymal <- colMeans(plot_data_mat[genes_mesenchymal_score,])
## make mesenchymal score
# genes_epithelial_score <- c("KRT19", "CLDN10", "GPX3", "SLC5A12")
genes_epithelial_score <- head(x = genes_plot_epithelial, n = 5); genes_epithelial_score
scores_epithelial <- colMeans(plot_data_mat[genes_epithelial_score,])

# make data frame ---------------------------------------------------------
emtscores_df <- data.frame(id_bycluster_bycellgroup_byaliquot = colnames(plot_data_mat),
                           Score_mesenchymal = scores_mesenchymal,
                           Score_epithelial = scores_epithelial)
emtscores_df <- emtscores_df %>%
  mutate(Sample_id = str_split_fixed(string = id_bycluster_bycellgroup_byaliquot, pattern = "_", n = 3)[,1]) %>%
  mutate(Cluster_id = str_split_fixed(string = id_bycluster_bycellgroup_byaliquot, pattern = "_", n = 3)[,2]) %>%
  select(-id_bycluster_bycellgroup_byaliquot) %>%
  select(Sample_id, Cluster_id, Score_mesenchymal, Score_epithelial) %>%
  mutate(EMT_potential = ifelse((scores_mesenchymal >= quantile(x = scores_mesenchymal, probs = 0.9)) & (scores_epithelial <= quantile(x = scores_epithelial, probs = 0.1)),
                                "High", "Low")) %>%
  arrange(EMT_potential, desc(Score_mesenchymal))

# write outputs -----------------------------------------------------------
file2write <- paste0(dir_out, "EMT_scores_by_seurattumorcluster.", run_id, ".tsv")
write.table(x = emtscores_df, file = file2write, sep = "\t", row.names = F, quote = F)

