# Yige Wu @WashU Feb 2020
## for plotting the genomics landscape of integrated single cell data

# source ------------------------------------------------------------------
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")

# set run id  ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)

# set output directory ----------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input seurat object master list -----------------------------------------
seurat_summary <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/summary/ccRCC_snRNA_Downstream_Processing - Seurat_Preprocessing.20200207.v1.tsv", data.table = F)
seurat_summary2process <- seurat_summary %>%
  filter(Cellranger_reference_type == "pre-mRNA") %>%
  filter(Proceed_for_downstream == "Yes") %>%
  filter(!(Aliquot %in% c("CPT0001540013", "CPT0002270013", "CPT0015810004", "CPT0023690004", "CPT0025110004", "CPT0063630004", "CPT0065690004", "CPT0000870003", "CPT0075720013"))) %>%
  mutate(Path_seurat_object = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                     "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                     "/", Aliquot, FACS, "_processed.rds"))
seurat_summary2process$Path_seurat_object


# input cluster to cell type table ----------------------------------------
cluster2celltype <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - Individual.AllCluster2Cell_Type.20200207.v2.tsv", data.table = F)

# Loop: for each aliquot, input seurat object and infercnv output, plot important genes on UMAP ---------
infercnv_run_id <- "20200207.v1"
dir_infercnv_all_runs <- "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/InferCNV/outputs/"
dir_infercnv_output <- paste0(dir_infercnv_all_runs, "Individual.", infercnv_run_id, "/")

# input infercnv results ------------------------------------------------
sup_tab_infercnv <- NULL
for (snRNA_aliquot_id_tmp in snRNA_aliquot_ids) {
  ## input infercnv outputs
  obs_cnv_state_mat <- fread(input = paste0(dir_infercnv_output, snRNA_aliquot_id_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"), data.table = F)
  ref_cnv_state_mat <- fread(input = paste0(dir_infercnv_output, snRNA_aliquot_id_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.references.txt"), data.table = F)
  
  # number_tumor_cells <- ncol(obs_cnv_state_mat) - 1
  # number_tumor_cells
  # 
  # number_normal_cells <- ncol(ref_cnv_state_mat) - 1
  # number_normal_cells
  
  ## filter the outputs to include only the CNV target genes
  obs_cnv_state_mat <- obs_cnv_state_mat %>%
    rename(gene_symbol = V1) %>%
    filter(gene_symbol %in% ccrcc_cna_genes_df$gene_symbol)
  
  ref_cnv_state_mat <- ref_cnv_state_mat %>%
    rename(gene_symbol = V1) %>%
    filter(gene_symbol %in% ccrcc_cna_genes_df$gene_symbol)
  
  ## transform wide data frame to long data frame 
  obs_cnv_state_mat.m <- melt(obs_cnv_state_mat, id.vars = c("gene_symbol"))
  ref_cnv_state_mat.m <- melt(ref_cnv_state_mat, id.vars = c("gene_symbol"))
  
  ## combine the CNV result for observation cell group and reference cell group
  cnv_state_mat.m <- rbind(obs_cnv_state_mat.m, ref_cnv_state_mat.m)
  
  ## filter out copy neutral results to save space
  cnv_state_mat.m <- cnv_state_mat.m %>%
    rename(cna_state = value) %>%
    filter(cna_state != 1)
  
  ## map barcode to cluster
  ### input individually processed seurat object
  seurat_obj_path <- seurat_summary2process$Path_seurat_object[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
  seurat_obj_path
  seurat_object <- readRDS(file = seurat_obj_path)
  ### get meta data
  meta_data_tmp <- seurat_object@meta.data
  meta_data_tmp$barcode <- rownames(meta_data_tmp)
  ## map barcode using the meta data
  cnv_state_mat.m$cluster <- meta_data_tmp[as.vector(cnv_state_mat.m$variable),"seurat_clusters"]
  cnv_state_mat.m$cluster <- as.numeric(as.vector(cnv_state_mat.m$cluster))
  rm(seurat_object)
  
  ## filter barcodes by only nephron epithelium
  ### get the cell type assignments for only this aliquot
  cluster2celltype_aliquot <- cluster2celltype %>%
    filter(Aliquot == snRNA_aliquot_id_tmp)
  ### get clusters that are assigned as nephron epithelium
  nep_epi_clusters <- cluster2celltype_aliquot$Cluster[cluster2celltype_aliquot$Most_Enriched_Cell_Group == "Nephron_Epithelium"]
  ### filter CNV result by cluster
  cnv_state_mat.m <- cnv_state_mat.m %>%
    filter(cluster %in% nep_epi_clusters)
  
  ## annotate the gene to expected copy state
  cnv_state_mat.m$gene_cna_type <- mapvalues(x = cnv_state_mat.m$gene_symbol, from = ccrcc_cna_genes_df$gene_symbol, to = as.vector(ccrcc_cna_genes_df$gene_cna_type))
  
  ## filter CNV results to only expected ones (for example, for 3p genes genes only keep deletion events)
  cnv_state_mat.m <- cnv_state_mat.m %>%
    filter((gene_cna_type == "Loss" & cna_state %in% c(0.5, 0)) | (gene_cna_type == "Gain" & cna_state %in% c(1.5, 2)))
  
  ## count the CNV events by cluster
  cnv_state_tab <- cnv_state_mat.m %>%
    select(gene_symbol, cluster) %>%
    table() %>%
    as.data.frame() %>%
    mutate(aliquot = snRNA_aliquot_id_tmp)
  
  ## combine with super table
  sup_tab_infercnv <- rbind(cnv_state_mat.m, sup_tab_infercnv)
}
sup_tab_infercnv <- sup_tab_infercnv %>%
  mutate(cna_text = paste0(gene_symbol, "_", gene_cna_type)) %>%
  mutate(sn_cna_cat = map_infercnv_state2category(copy_state = cna_state))

# load meta data ----------------------------------------------------------
meta_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/sample_info/make_meta_data/meta_data.20190924.v1.tsv", data.table = F)
sup_tab <- merge(sup_tab_infercnv, meta_tab, by.x = c("aliquot"), by.y = c("Specimen.ID.snRNA"), all.x = T)
sup_tab$Case.ID[sup_tab$aliquot == "CPT0075140002"] <- "C3N-01200"
sup_tab$Case.ID[sup_tab$aliquot == "CPT0025890002"] <- "C3N-00733"


# add chromosomal region --------------------------------------------------
sup_tab$chr_region <- mapvalues(x = sup_tab$gene_symbol, from = ccrcc_cna_genes_df$gene_symbol, to = as.vector(ccrcc_cna_genes_df$chr_region))
sup_tab$chr_region <- factor(sup_tab$chr_region, levels = c("3p", "5q", "14q",
                                                            "9p21", 
                                                            "10q23",
                                                            "1p31",
                                                            "6q24",
                                                            "3p12",
                                                            "9p23",
                                                            "14q24",
                                                            "3p26",
                                                            "1q32",
                                                            "8q24",
                                                            "9q24"))

# plot heatmap showing sn-based CNA Frequency------------------------------------------------------------
tab2p <- sup_tab
tab2p$perc_cna_in_all_cell[tab2p$perc_cna_in_all_cell == 0] <- NA

p <- ggplot()
p <- p + geom_point(data = tab2p, mapping = aes(x = gene_symbol, y = aliquot, size = perc_cna_in_all_cell, fill = sn_cna_cat), shape = 21)
p <- p + scale_fill_manual(values = copy_number_colors)
p <- p + facet_grid(Case.ID~chr_region + gene_cna_type, scales = "free", space = "free", shrink = T)
p <- p + theme_bw()
p <- p + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
p <- p + theme(panel.spacing = unit(0, "lines"))
p
png(file = paste0(dir_out, "TCGA_CNA_Genes_snCNA.withlegend", ".", run_id, ".png"), width = 1200, height = 800, res = 150)
print(p)
dev.off()

# write out CNA frequency table -------------------------------------------
write.table(x = tab2p, file = paste0(dir_out, "TCGA_CNA_Genes_snCNA_Frequency.tsv"), quote = F, sep = "\t", row.names = F)
