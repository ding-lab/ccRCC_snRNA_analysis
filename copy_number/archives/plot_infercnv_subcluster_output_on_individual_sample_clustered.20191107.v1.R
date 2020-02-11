# Yige Wu @WashU Oct 2019

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")

# set parameters ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
integration_id <- "20191021.v1"
dir_infercnv_output <- "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/InferCNV/outputs/"

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

# input cluster2celltype table --------------------------------------------
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - Individual.AllCluster2Cell_Type.20191106.v1.tsv", data.table = F)

# input the HIF downstream to plot ----------------------------------------
HIF_targets2plot <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/proteomics/plot_HIF_pahtway_protein_heatmap/ccRCC_snRNA_Downstream_Processing - HIF_Target_Summary.tsv", data.table = F)

# set genes to plot -------------------------------------------------------
genes2plot <- ccrcc_cna_genes_df$gene_symbol
genes2plot <- "VHL"
genes2plot <- "HIF1A"
genes2plot <- c("VHL", "HIF1A", "EPAS1")
genes2plot <- c("PBRM1", "BAP1", "SETD2")
genes2plot <- c("VHL", "HIF1A", "EPAS1", "PBRM1", "BAP1", "SETD2")
genes2plot <- c("ARNT", "EPAS1")

# plot by each aliquot ----------------------------------------------------
# snRNA_aliquot_id_tmp <- c("CPT0019130004")
# snRNA_aliquot_id_tmp <- c("CPT0025890002")
# snRNA_aliquot_id_tmp <- c("CPT0086350004")
# snRNA_aliquot_id_tmp <- c("CPT0010110013")
# snRNA_aliquot_id_tmp <- c("CPT0001180011")
# snRNA_aliquot_id_tmp <- c("CPT0001260013")
# snRNA_aliquot_id_tmp <- c("CPT0075140002")
snRNA_aliquot_id_tmp <- c("CPT0001220012")
snRNA_aliquot_id_tmp <- "CPT0001220012"

for (snRNA_aliquot_id_tmp in snRNA_aliquot_ids) {
  # input individually processed seurat object ---------------------------------------------
  seurat_obj_path <- seurat_summary2process$Path_seurat_object[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
  seurat_obj_path
  seurat_object <- readRDS(file = seurat_obj_path)
  
  # get umap coordinates ----------------------------------------------------
  umap_tab <- FetchData(seurat_object, vars = c("orig.ident", "ident", "UMAP_1", "UMAP_2"))
  umap_tab$barcode <- rownames(umap_tab)
  
  ## get the cluster2celltype info for current aliquot
  cluster2celltype_tab_tmp <- cluster2celltype_tab %>%
    filter(Aliquot == snRNA_aliquot_id_tmp)
  seurat_obj@meta.data$cluster_cell_type <- plyr::mapvalues(seurat_obj@meta.data$seurat_clusters, from = cluster2celltype_tab_tmp$Cluster, to = cluster2celltype_tab_tmp$Enriched_Cell_Type_Abbr)
  
  p <- DimPlot(seurat_obj, reduction = "umap", group.by = "cluster_cell_type", label = T, label.size	= 5, repel = T)
  label_data <- p$layers[[2]]$data
  label_data <- label_data %>%
    filter(cluster_cell_type %in% cluster2celltype_tab_tmp$Enriched_Cell_Type_Abbr[cluster2celltype_tab_tmp$Is_Malignant == "Yes"])
  
  # input infercnv observations ---------------------------------------------
  tumor_cnv_state_mat <- fread(input = paste0(dir_infercnv_output, "integration.", integration_id, "/", snRNA_aliquot_id_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.observations.txt"), data.table = F)
  ref_cnv_state_mat <- fread(input = paste0(dir_infercnv_output, "integration.", integration_id, "/", snRNA_aliquot_id_tmp, "/infercnv.14_HMM_predHMMi6.rand_trees.hmm_mode-subclusters.Pnorm_0.5.repr_intensities.references.txt"), data.table = F)
  dim(tumor_cnv_state_mat)
  dim(ref_cnv_state_mat)
  
  cnv_state_df <- rbind(melt(tumor_cnv_state_mat, id.vars = c("V1")), melt(ref_cnv_state_mat, id.vars = c("V1")))
  
  for (gene_tmp in genes2plot) {
    if (gene_tmp %in% cnv_state_df$V1) {
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
      
      tab2p <- tab2p %>%
        arrange(desc(cnv_cat))
      
      
      # plot a UMAP plot for copy number metric ---------------------------------
      p <- ggplot() +
        geom_point(data = tab2p, mapping = aes(UMAP_1, UMAP_2, color=cnv_cat), alpha = 1, size = 0.3) +
        scale_color_manual(values = copy_number_colors) +
        ggtitle(paste0(snRNA_aliquot_id_tmp, "_", gene_tmp, "_Copy_Number_Status"))
      p <- p + geom_text_repel(data = label_data, mapping = aes(UMAP_1, UMAP_2, label = cluster_cell_type))
      p <- p + theme_bw() +
        theme(panel.border = element_blank(), panel.grid.major = element_blank(),
              panel.grid.minor = element_blank())
      p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                              axis.text.y=element_blank(),axis.ticks=element_blank(),
                              axis.title.x=element_blank(),
                              axis.title.y=element_blank())
      file2write <- paste(dir_out, snRNA_aliquot_id_tmp,"_Individual_Clustered_", gene_tmp, ".FeaturePlot_CNA.", run_id, ".png", sep="")
      png(file2write, width = 1100, height = 800, res = 150)
      print(p)
      dev.off()
      
      
    } else {
      print(paste0(gene_tmp, " not in the infercnv result!"))
    }
  }
}

# plot a Dotplot for copy number metric ---------------------------------
genes2plot <- HIF_targets2plot$Gene_Symbol[HIF_targets2plot$Tumor == 1 & is.na(HIF_targets2plot$Stromal) & is.na(HIF_targets2plot$Immune)]
seurat_obj@meta.data$VHL_cna_cat <- plyr::mapvalues(rownames(seurat_obj@meta.data), from = tab2p$barcode, to = tab2p$cnv_cat)
seurat_obj@meta.data$cluster_cell_type <- plyr::mapvalues(seurat_obj@meta.data$seurat_clusters, from = cluster2celltype_tab_tmp$Cluster, to = cluster2celltype_tab_tmp$Enriched_Cell_Type_Abbr)
seurat_obj@meta.data$cluster_is_malignant <- plyr::mapvalues(seurat_obj@meta.data$seurat_clusters, from = cluster2celltype_tab_tmp$Cluster, to = cluster2celltype_tab_tmp$Is_Malignant)

seurat_obj2 <- subset(seurat_obj, subset = (VHL_cna_cat %in% c("Neutral", "Loss of one copy")))

p <- DotPlot(seurat_obj2, features = genes2plot, split.by = "VHL_cna_cat") + RotatedAxis()
p$data$Cluster <- str_split_fixed(string = p$data$id, pattern = "_", n = 2)[,1]
p$data$cluster_is_malignant <- ifelse(p$data$Cluster %in% cluster2celltype_tab_tmp$Cluster[cluster2celltype_tab_tmp$Is_Malignant == "Yes"], "Malignant_Cell_Clusters", "Nom-malignant_Cell_Clusters")
p <- p + facet_grid(cluster_is_malignant~., scales = "free", space = "free", drop = T)
file2write <- paste(dir_out, snRNA_aliquot_id_tmp,"_Individual_Clustered_VHL_CNA_vs_HIF_Targets", ".Dotplot.", run_id, ".png", sep="")
png(file2write, width = 1200, height = 800, res = 150)
print(p)
dev.off()

# plot a Dotplot for copy number metric ---------------------------------
p <- VlnPlot(seurat_obj2, features = "VEGFA", split.by = "VHL_cna_cat")
p$layers[[2]]$aes_params$size <- 0.5
p$layers[[2]]$aes_params$alpha <- 0.3
p <- p + scale_x_discrete("Cluster", labels = mapvalues(x = ggplot_build(p)$layout$panel_params[[1]]$x.labels, from = cluster2celltype_tab_tmp$Cluster, to = cluster2celltype_tab_tmp$Enriched_Cell_Type_Abbr))
file2write <- paste(dir_out, snRNA_aliquot_id_tmp,"_Individual_Clustered_VHL_CNA_vs_HIF_Targets", ".Vlnplot.", run_id, ".png", sep="")
png(file2write, width = 1200, height = 1000, res = 150)
print(p)
dev.off()



