# Yige Wu @ WashU 20189 Aug
## write the tables for cell type classification project


# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")

# Set version -----------------------------------------------------------------
version_num <- 1


# input seurat onbjects ---------------------------------------------------
qced_object <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/scRNA/C3N-01200-02_normalref/C3N-01200-02_normalref_Clustered_seurat_object_nFeature_over200_under5000_nCount_over1000_under50000_percent.mito_0.1.RDS")
sample_id <- "CPT0075130004"

# write the meta_data table  -----------------------------------------------------
meta_data_tab <- qced_object@meta.data
meta_data_tab$cell_barcode <- rownames(meta_data_tab)
meta_data_tab$sample_id <- sample_id
meta_data_tab$technique <- "sn"

cell_type_manual <- data.frame(seurat_clusters = c(0, 
                                                   1, 
                                                   2, 
                                                   3, 
                                                   4, 
                                                   5, 
                                                   6,
                                                   7,
                                                   8,
                                                   9,
                                                   10),
                               cell_type = c("Unknown", 
                                             "ClearCellRenalCellCarcinoma", 
                                             "MacroMonoNeutro", 
                                             "ClearCellRenalCellCarcinoma", 
                                             "Tcell", 
                                             "Unknown", 
                                             "Unknown",
                                             "Unknown",
                                             "Unknown",
                                             "Endothelial",
                                             "Unknown"))
meta_data_tab <- merge(meta_data_tab, cell_type_manual, by = c("seurat_clusters"), all.x = T)
meta_data_tab %>%
  select(cell_type) %>%
  table()
meta_data_tab4pollock <- meta_data_tab %>%
  select(sample_id, cell_barcode, cell_type, technique)
write.table(x = meta_data_tab4pollock, file = paste0(makeOutDir(), sample_id, "_meta_data", ".", format(Sys.Date(), "%Y%m%d") , ".v", version_num, ".txt"), quote = F, sep = "\t", row.names = F)

# write the raw expression table ------------------------------------------
raw_expression_tab <-  as.data.frame(qced_object@assays$RNA@counts)
raw_expression_tab %>%
  head()
rownames(raw_expression_tab) %>%
  head()
write.table(x = raw_expression_tab, file = paste0(makeOutDir(), sample_id, "_raw_expression_tab", ".", format(Sys.Date(), "%Y%m%d") , ".v", version_num, ".txt"), quote = F, sep = "\t", row.names = F)


