# Yige Wu @WashU Nov 2019
## for plotting the marker genes for integrated object, showing cell of origin

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

# input marker table ------------------------------------------------------
gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Kidney_Markers/RCC_marker_gene_table_and_literature_review - Gene2CellType_Tab_w.HumanProteinAtlas.20191024.v1.tsv")


# input cluster 2 cell type table -----------------------------------------
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - Individual.AllCluster2Cell_Type.20191112.v1.tsv", data.table = F)

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


# plot dotplot by sample --------------------------------------------------
snRNA_aliquot_id_tmp <- "CPT0001220012"

cluster2celltype_w_mut_tab <- NULL
for (snRNA_aliquot_id_tmp in snRNA_aliquot_ids) {
  ## get the cluser 2 cell table
  cluster2celltype_tab_tmp <- cluster2celltype_tab %>%
    filter(Aliquot == snRNA_aliquot_id_tmp)
  
  ## input 10XMapping results
  mutation_map_tab <- fread(input = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/10Xmapping/outputs/", snRNA_aliquot_id_tmp, "/", snRNA_aliquot_id_tmp, "_mapping_heatmap_0.txt"), data.table = F)
  
  ## get variants mapped in 10XMapping results
  mutation_map_tab <- mutation_map_tab %>%
    mutate(gene_symbol = str_split_fixed(string = Mutatation, pattern = "-", n = 3)[,1])
  mutation_map_tab.m <- melt(mutation_map_tab, id.vars = c("gene_symbol", "Mutatation"))
  mutation_map_tab.m <- mutation_map_tab.m %>%
    mutate(allele_type = str_split_fixed(string = Mutatation, pattern = "-", n = 3)[,3])
  
  mutation_map_tab.var <- mutation_map_tab.m %>%
    filter(allele_type == "Var") %>%
    filter(!is.na(value) & value > 0) %>%
    rename(barcode = variable) %>%
    mutate(mutation = paste0(gene_symbol, "-", str_split_fixed(string = Mutatation, pattern = "-", n = 3)[,2]))
  
  ## input seurat object
  seurat_obj_path <- seurat_summary2process$Path_seurat_object[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
  seurat_obj <- readRDS(file = seurat_obj_path)
  
  ## get meta data table
  meta_data_tab_tmp <- seurat_obj@meta.data
  meta_data_tab_tmp$barcode <- rownames(meta_data_tab_tmp)
  
  ## map barcode to cluster in mutation mapping table
  mutation_map_tab.var$cluster <- mapvalues(x = mutation_map_tab.var$barcode, from = meta_data_tab_tmp$barcode, to = as.vector(meta_data_tab_tmp$seurat_clusters))
  
  ## summarize all the mutation mapped to each cluster
  cluster2celltype_tab_tmp$Somatic_Mutation_Mapped <- sapply(X = cluster2celltype_tab_tmp$Cluster, FUN = function(c, mut_tab) {
    mutations_tmp <- unique(mut_tab$mutation[mut_tab$cluster == c])
    mutation_text_tmp <- paste0(mutations_tmp, collapse = ";")
    return(mutation_text_tmp)
  }, mut_tab = mutation_map_tab.var)
  
  ## merge mutaiton mapping with super table
  cluster2celltype_w_mut_tab <- rbind(cluster2celltype_tab_tmp, cluster2celltype_w_mut_tab)
}

write.table(x = cluster2celltype_w_mut_tab, file = paste0(dir_out, "Individual.AllCluster2Cell_Type.MutationMapped.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)
