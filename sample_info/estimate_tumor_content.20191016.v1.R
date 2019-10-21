# Yige Wu @WashU Oct 2019
## for plotting the fraction of immune cell populations across samples

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
object2plot <- readRDS("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/analysis_results/integration/integrate_seurat_objects/20191015.v1/Renal_Integrated.20191015.v1.RDS")
DefaultAssay(object2plot) <- "RNA"
aliquot_ids <- unique(object2plot@meta.data$orig.ident)

# input cluster cell type assignment --------------------------------------
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - AllCluster2Cell_Type.20191016.v1.tsv", data.table = F)

# get number of barcodes per cluster --------------------------------------
sn_cell_num_tab <- data.frame(object2plot@meta.data %>%
                                select(orig.ident, seurat_clusters) %>%
                                table())

colnames(sn_cell_num_tab) <- c("snRNA_Aliquot_ID", "Cluster", "Num_Cluster_Barcode")
sn_cell_num_tab <- merge(sn_cell_num_tab, cluster2celltype_tab, by = c("Cluster"), all.x = T)

sn_cell_sum_tab <- sn_cell_num_tab %>%
  group_by(snRNA_Aliquot_ID) %>%
  summarise(Num_Aliquot_All_Cells = sum(Num_Cluster_Barcode))

sn_tumor_cell_sum_tab <- sn_cell_num_tab %>%
  filter(Is_Malignant == "Yes") %>%
  group_by(snRNA_Aliquot_ID) %>%
  summarise(Num_Aliquot_Tumor_Cells = sum(Num_Cluster_Barcode))

sn_tumor_cell_sum_tab <- merge(sn_tumor_cell_sum_tab, sn_cell_sum_tab, by = c("snRNA_Aliquot_ID"), all.x = T)
sn_tumor_cell_sum_tab <- sn_tumor_cell_sum_tab %>%
  mutate(Perc_Aliquot_Tumor_Cell = Num_Aliquot_Tumor_Cells/Num_Aliquot_All_Cells)

# input meta data ---------------------------------------------------------
meta_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/analysis_results/sample_info/make_meta_data/meta_data.20190924.v1.tsv", data.table = F)
meta_tab_in_use <- meta_tab %>%
  dplyr::filter(meta_tab$Specimen.ID.snRNA %in% aliquot_ids)
rownames(meta_tab_in_use) <- meta_tab_in_use$Case.ID


# input tumor purity ESTIMATE-RNA based result ---------------------------------------------------
estimate_tab <- readxl::read_excel("./Ding_Lab/Projects_Current/CPTAC/PGDAC/ccRCC_discovery_manuscript/ccRCC Manuscript/CPTAC3-ccRCC-SupplementaryTables_Final/Table S7.xlsx", sheet = "ESTIMATE scores")
estimate_purity_rna <- as.numeric(as.vector(as.data.frame(estimate_tab[7,2:176])))
estimate_purity_rna_df <- data.frame(ESTIMATE_TumorPurity_RNA = estimate_purity_rna, Specimen.ID.bulk = unlist(estimate_tab[3,2:176]))
estimate_purity_rna_df <- estimate_purity_rna_df %>%
  filter(Specimen.ID.bulk %in% meta_tab_in_use$Specimen.ID.bulk)
estimate_purity_rna_df$Specimen.ID.snRNA <- mapvalues(estimate_purity_rna_df$Specimen.ID.bulk, from = meta_tab_in_use$Specimen.ID.bulk, to = meta_tab_in_use$Specimen.ID.snRNA)

# merge ESTIMATE tumor content into sn-based tumor content ------------------------------
sn_tumor_cell_sum_tab$ESTIMATE_TumorPurity_RNA <- mapvalues(sn_tumor_cell_sum_tab$snRNA_Aliquot_ID, from = estimate_purity_rna_df$Specimen.ID.snRNA, to = estimate_purity_rna_df$ESTIMATE_TumorPurity_RNA)

file2write <- paste0(dir_out, "Perc_Tumor_Content_w_ESTIMATE.", run_id, ".tsv")
write.table(x = sn_tumor_cell_sum_tab, file = file2write, quote = F, row.names = F, sep = "\t")




