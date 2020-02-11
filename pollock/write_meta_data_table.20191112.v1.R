# Yige Wu @WashU Nov 2019
## for making the meta data table for pollock training of ccRCC snRNA-seq data

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


# set aliquot ids to be processed -----------------------------------------
snRNA_aliquot_ids <- c("CPT0019130004", "CPT0001260013", "CPT0086350004", "CPT0010110013", "CPT0001180011", "CPT0025890002", "CPT0075140002", "CPT0020120013", "CPT0001220012", "CPT0014450005")

# input seurat processing summary ------------------------------------------------
seurat_summary <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/summary/ccRCC_snRNA_Downstream_Processing - Seurat_Preprocessing.20191021.v1.tsv", data.table = F)
seurat_summary2process <- seurat_summary %>%
  filter(Cellranger_reference_type == "pre-mRNA") %>%
  filter(Proceed_for_downstream == "Yes") %>%
  filter(Aliquot %in% snRNA_aliquot_ids) %>%
  mutate(Path_seurat_object = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                     "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                     "/", Aliquot, FACS, "_processed.rds")) %>%
  mutate(Paht_deg_table = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                 "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                 "/", Aliquot, FACS, ".DEGs.Pos.txt"))
seurat_summary2process$Path_seurat_object

# input cluster 2 cell type table -----------------------------------------
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - Individual.AllCluster2Cell_Type.20191112.v1.tsv", data.table = F)

# make meta data table --------------------------------------------------
## sample_id, cell_barcode, cell_type, technique(sc or sn)
meta_data_tab <- NULL

# snRNA_aliquot_id_tmp <- "CPT0001180011"
for (snRNA_aliquot_id_tmp in snRNA_aliquot_ids) {
  ## input seurat object
  seurat_obj_path <- seurat_summary2process$Path_seurat_object[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
  seurat_obj <- readRDS(file = seurat_obj_path)
  
  ## get cluster 2 cell type table
  cluster2celltype_tab_tmp <- cluster2celltype_tab %>%
    filter(Aliquot == snRNA_aliquot_id_tmp)
  
  ## make meta data table
  meta_data_tab_tmp <- data.frame(sample_id = snRNA_aliquot_id_tmp, 
                                  cell_barcode = rownames(seurat_obj@meta.data),
                                  cluster = seurat_obj@meta.data$seurat_clusters,
                                  technique = "sn")
  
  meta_data_tab_tmp$cell_type <- mapvalues(x = meta_data_tab_tmp$cluster, from = cluster2celltype_tab_tmp$Cluster, to = cluster2celltype_tab_tmp$Enriched_Cell_Type_Abbr)
  
  meta_data_tab <- rbind(meta_data_tab_tmp %>%
                           select(sample_id, cell_barcode, cell_type,technique), meta_data_tab)
}
write.table(x = meta_data_tab, file = paste0(dir_out, "ccRCC_meta_data.", run_id, ".tsv"), quote = F, row.names = F, sep = "\t")

