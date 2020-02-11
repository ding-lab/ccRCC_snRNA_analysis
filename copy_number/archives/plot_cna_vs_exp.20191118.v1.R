# Yige Wu @WashU Nov 2019
## make dotplot and violin plot for showing gene expression in CNV group vs CN neutral group

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


# set the integration id to locate the infercnv run -----------------------
integration_id <- "20191021.v1"


# set aliquot id ----------------------------------------------------------
snRNA_aliquot_ids <- c("CPT0019130004", "CPT0001260013", "CPT0086350004", "CPT0010110013", "CPT0001180011", "CPT0025890002", "CPT0075140002", "CPT0020120013", "CPT0001220012", "CPT0014450005")

# input seurat processing summary ------------------------------------------------
seurat_summary <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/summary/ccRCC_snRNA_Downstream_Processing - Seurat_Preprocessing.20191021.v1.tsv", data.table = F)
seurat_summary2process <- seurat_summary %>%
  filter(Cellranger_reference_type == "pre-mRNA") %>%
  filter(Proceed_for_downstream == "Yes") %>%
  filter(Aliquot %in% snRNA_aliquot_ids) %>%
  mutate(Path_seurat_object = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                     "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                     "/", Aliquot, FACS, "_processed.rds"))
seurat_summary2process$Path_seurat_object


# input cluster 2 cell type table -----------------------------------------
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - Individual.AllCluster2Cell_Type.20191112.v1.tsv", data.table = F)

# set genes to plot -------------------------------------------------------
genes2plot <- ccrcc_cna_genes_df$gene_symbol

# plot by each aliquot ----------------------------------------------------
snRNA_aliquot_id_tmp <- "CPT0001220012"

for (snRNA_aliquot_id_tmp in snRNA_aliquot_ids) {
  # input individually processed seurat object ---------------------------------------------
  seurat_obj_path <- "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/individual_cluster/recluster_tumor/20191109.v1/CPT0001220012.Malignant_Reclustered.20191109.v1.RDS"
  seurat_obj_path
  seurat_object <- readRDS(file = seurat_obj_path)
  Idents(seurat_object) <- rownames(seurat_object@meta.data)
  
  # input infercnv observations ---------------------------------------------
  tumor_cnv_state_mat <- fread(input = paste0(dir_infercnv_output, "integration.", integration_id, "/", snRNA_aliquot_id_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"), data.table = F)
  ref_cnv_state_mat <- fread(input = paste0(dir_infercnv_output, "integration.", integration_id, "/", snRNA_aliquot_id_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.references.txt"), data.table = F)
  dim(tumor_cnv_state_mat)
  dim(ref_cnv_state_mat)
  
  cnv_state_df <- rbind(melt(tumor_cnv_state_mat, id.vars = c("V1")), melt(ref_cnv_state_mat, id.vars = c("V1")))

  ## for each gene, divide the tumor cells into 2 groups, one with CNA, one without
  for (gene_tmp in genes2plot) {
    if (gene_tmp %in% cnv_state_df$V1) {
      infercnv_observe_gene_tab <- cnv_state_df %>%
        rename(gene_symbol = V1) %>%
        filter(gene_symbol == gene_tmp) %>%
        mutate(barcode = str_split_fixed(string = variable, pattern = "_", n = 2)[,1]) %>%
        rename(copy_state = value) %>%
        select(gene_symbol, barcode, copy_state)
      
      infercnv_observe_gene_tab$cna_group <- ifelse(infercnv_observe_gene_tab$copy_state == 1.0, "Neutral",
                                                    ifelse(infercnv_observe_gene_tab$copy_state < 1.0, "Loss", "Gain"))
      
      object2plot <- seurat_object
      object2plot@meta.data$cna_group <- mapvalues(x = rownames(object2plot@meta.data), from = infercnv_observe_gene_tab$barcode, to = as.vector(infercnv_observe_gene_tab$cna_group))
      malignant_clusters <- cluster2celltype_tab$Cluster[cluster2celltype_tab$Aliquot == snRNA_aliquot_id_tmp]
      object2plot@meta.data$is_malignant <- ifelse(object2plot@meta.data$seurat_clusters %in% malignant_clusters, T, F)
      
      malignant_gene_cna <- ccrcc_cna_genes_df[ccrcc_cna_genes_df$gene_symbol == gene_tmp, "gene_cna_type"] %>% as.vector()
      
      object2plot <- subset(object2plot, subset = (is_malignant == T & (cna_group %in% c("Neutral", malignant_gene_cna))))
      
      object2plot@meta.data$Cluster_cna_cat <- paste0("Cluster", object2plot@meta.data$seurat_clusters, "\n", gene_tmp, "_", object2plot@meta.data$cna_group)
      
      p <- VlnPlot(object2plot, features = gene_tmp, group.by = "Cluster_cna_cat", split.by = "cna_group")
      
      
      p <- DotPlot(object2plot, features = gene_tmp, group.by = "Cluster_cna_cat", col.max = 2.5, col.min = 0) + RotatedAxis()
      p$data$Cluster <- str_split_fixed(string = p$data$id, pattern = "_|\\n", n = 2)[,1]
      p <- p + facet_grid(Cluster~., scales = "free", space = "free", drop = T)
      p <- p + theme(panel.spacing = unit(0, "lines"),
                     strip.background = element_blank(),
                     panel.border = element_rect(colour = "black"),
                     strip.text.y = element_text(vjust = 0.5, size = 12, face = "bold"),
                     strip.placement = "outside")
      file2write <- paste(dir_out, snRNA_aliquot_id_tmp,".Individual_Malignant_Reclustered.", gene_tmp, ".Loss.", "VEGFA", ".Dotplot.", run_id, ".png", sep="")
      png(file2write, width = 1200, height = 1000, res = 150)
      print(p)
      dev.off()
      
      stop("m")
      
    } else {
      print(paste0(gene_tmp, " not in the infercnv result!"))
    }
  }
}


# plot a Violin plot for copy number metric ---------------------------------
gene_tmp <- "VHL"
genes2plot <- HIF_targets2plot$Gene_Symbol[HIF_targets2plot$Tumor == 1 & is.na(HIF_targets2plot$Stromal) & is.na(HIF_targets2plot$Immune)]

for (gene_tmp in c("PBRM1")) {
  infercnv_observe_gene_tab <- cnv_state_df %>%
    rename(gene_symbol = V1) %>%
    filter(gene_symbol == gene_tmp) %>%
    mutate(barcode = str_split_fixed(string = variable, pattern = "_", n = 2)[,1]) %>%
    rename(copy_state = value) %>%
    select(gene_symbol, barcode, copy_state)
  tab2p <- umap_tab
  tab2p <- merge(tab2p, infercnv_observe_gene_tab, by = c("barcode"), all.x = T)
  tab2p$cnv_cat <- map_infercnv_state2category(copy_state = tab2p$copy_state)
  tab2p$cnv_cat %>% table()
  
  seurat_object@meta.data$cna_cat <- plyr::mapvalues(rownames(seurat_object@meta.data), from = tab2p$barcode, to = tab2p$cnv_cat)
  seurat_object@meta.data$Cluster_cna_cat <- paste0("Cluster", seurat_object@meta.data$seurat_clusters, "\n", gene_tmp, "_", seurat_object@meta.data$cna_cat)
  
  seurat_obj2 <- subset(seurat_object, subset = (cna_cat %in% c("Neutral", "Loss of one copy")))
  
  p <- VlnPlot(seurat_obj2, features = "VEGFA", group.by = "Cluster_cna_cat", split.by = "cna_cat")
  p$layers[[1]]$geom_params$draw_quantiles <- T
  p$layers[[2]]$aes_params$size <- 0.5
  p$layers[[2]]$aes_params$alpha <- 0.3
  p$theme$axis.text.x$angle <- 90
  file2write <- paste(dir_out, snRNA_aliquot_id_tmp,".Individual_Malignant_Reclustered.", gene_tmp, ".Loss.", "VEGFA", ".Vlnplot.", run_id, ".png", sep="")
  png(file2write, width = 1500, height = 1000, res = 150)
  print(p)
  dev.off()
}

