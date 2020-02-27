# Yige Wu @WashU Feb 2020
## for each individual sample, isolating cells within the immune cell clusters and re-do clustering

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# input dependencies ------------------------------------------------------
## input seurat object master list
seurat_summary <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/summary/ccRCC_snRNA_Downstream_Processing - Seurat_Preprocessing.20200207.v1.tsv", data.table = F)
seurat_summary2process <- seurat_summary %>%
  filter(Cellranger_reference_type == "pre-mRNA") %>%
  filter(Proceed_for_downstream == "Yes") %>%
  mutate(Path_seurat_object = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                     "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                     "/", Aliquot, FACS, "_processed.rds"))
seurat_summary2process$Path_seurat_object
## input the cluster to cell type table
cluster2celltype_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - Individual.AllCluster2Cell_Type.20200220.v2.tsv", data.table = F)
## input gene to cell type info
gene2celltype_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Kidney_Markers/Gene2CellType_Tab.20200220.v1.tsv", data.table = F)


for (aliquot_tmp in seurat_summary2process$Aliquot) {
  # run reclustering by each aliquot ----------------------------------------
  ## create output directory by aliquot
  dir_out1 <- paste0(dir_out, aliquot_tmp, "/")
  dir.create(dir_out1)
  
  ## check if the reclustered object has been saved for this aliquot
  file2write <- paste0(dir_out1, aliquot_tmp, ".lymphoid_reclustered.", run_id, ".RDS")
  if (!file.exists(file2write)) {
    ## input individually processed seurat object
    seurat_obj_path <- seurat_summary2process$Path_seurat_object[seurat_summary2process$Aliquot == aliquot_tmp]
    seurat_obj_path
    seurat_object <- readRDS(file = seurat_obj_path)
    
    ## get the immune clusters for this aliquot
    clusters2process <- cluster2celltype_df$Cluster[cluster2celltype_df$Aliquot == aliquot_tmp & cluster2celltype_df$Most_Enriched_Cell_Type1 == "Lymphoid lineage immune cells"]
    clusters2process
    
    ## subset data
    new.object <- subset(seurat_object, idents = clusters2process)
    
    ## Run the standard workflow for clustering and visualization
    new.object <- FindVariableFeatures(object = new.object, selection.method = "vst", nfeatures = 2000)
    new.object <- ScaleData(new.object, verbose = F)
    new.object <- RunPCA(new.object, npcs = 30, verbose = FALSE)
    new.object <- RunUMAP(new.object, reduction = "pca", dims = 1:30)
    new.object <- FindNeighbors(new.object, reduction = "pca", dims = 1:30, force.recalc = T)
    new.object <- FindClusters(new.object, resolution = 0.5)
    saveRDS(object = new.object, file = file2write)
    
    
    # plot the cells back the original object dimensions ----------------------
    ### get the umap coordates for the original object and barcodes, and old cluster info
    plot_data_df <- FetchData(seurat_object, vars = c("ident", "UMAP_1", "UMAP_2"))
    plot_data_df$barcode <- rownames(plot_data_df)
    p <- DimPlot(seurat_object, reduction = "umap", label = T, label.size	= 5, repel = T)
    label_data_df <- p$layers[[2]]$data
    
    ### get the barcodes for new object, with new cluster info
    new_barcode_info_df <- new.object@meta.data
    new_barcode_info_df$barcode <- rownames(new_barcode_info_df)
    new_barcode_info_df <- new_barcode_info_df %>%
      select(barcode, seurat_clusters) %>%
      rename(new_cluster = seurat_clusters)
    plot_data_df <- merge(plot_data_df, new_barcode_info_df, by = c("barcode"), all.x = T)
    
    ## ggplot
    p <- ggplot()
    p <- p + geom_point(data = plot_data_df[is.na(plot_data_df$new_cluster),], 
                        mapping = aes(x = UMAP_1, y = UMAP_2), alpha = 0.3, size = 0.3, color = "grey50")
    p <- p + geom_point(data = plot_data_df[!is.na(plot_data_df$new_cluster),], 
                        mapping = aes(x = UMAP_1, y = UMAP_2, color=new_cluster), alpha = 0.8, size = 0.8)
    p <- p + geom_text_repel(data = label_data_df, mapping = aes(UMAP_1, UMAP_2, label = ident), size = 5, colour = "white")
    p <- p +
      theme_bw() +
      theme(panel.border = element_blank(), 
            panel.background = element_rect(fill = "black", colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
    p <- p + theme(legend.position = "right",
                   legend.text = element_text(size = 9))
    p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                            axis.text.y=element_blank(),axis.ticks=element_blank(),
                            axis.title.x=element_blank(),
                            axis.title.y=element_blank())
    ## save plot
    file2write <- paste(dir_out1, aliquot_tmp, ".lymphoid_reclustered.", "mapped_back.", run_id, ".png", sep="")
    png(file2write, width = 1000, height = 800, res = 150)
    print(p)
    dev.off()
    
    # dotplot -----------------------------------------------------------------
    if (length(unique(new.object@meta.data$seurat_clusters)) > 1) {
      genes2plot <-  intersect(gene2celltype_df$Gene, new.object@assays$RNA@counts@Dimnames[[1]])
      genes2plot <- unique(genes2plot)
      #### get the pct expressed for each gene in each cluster
      p <- DotPlot(object = new.object, features = genes2plot, col.min = 0)
      plot_data <- p$data
      #### transform the dataframe to matrix to better filter out genes with too low expressin
      plot_matrix <- dcast(data = plot_data, formula = features.plot ~ id, value.var = "pct.exp")
      plot_matrix %>% head()
      #### filter for genes that are expressed in >25% of one cluster at least
      genes2plot_filtered <- as.vector(plot_matrix[rowSums(plot_matrix[,unique(as.vector(plot_data$id))] > 10) >= 1, "features.plot"])
      ## dotplot
      p <- DotPlot(object = new.object, features = genes2plot_filtered, col.min = 0)
      p$data$gene_cell_type_group <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type_Group)
      p$data$gene_cell_type1 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type1)
      p$data$gene_cell_type2 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type2)
      p$data$gene_cell_type3 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type3)
      p$data$gene_cell_type4 <- plyr::mapvalues(p$data$features.plot, from = gene2celltype_df$Gene, to = gene2celltype_df$Cell_Type4)
      p <- p+ RotatedAxis()
      p <- p + facet_grid(.~gene_cell_type_group + gene_cell_type1 + gene_cell_type2 + gene_cell_type3 + gene_cell_type4, scales = "free", space = "free", drop = T)
      p <- p + theme(panel.spacing = unit(0, "lines"),
                     strip.background = element_blank(),
                     panel.border = element_rect(colour = "black"),
                     panel.grid.major = element_line(colour = "grey50"),
                     strip.text.x = element_text(angle = 0, vjust = 0.5),
                     axis.text.x = element_text(size = 15, face = "bold"),
                     strip.placement = "outside")
      
      ## save plot
      file2write <- paste0(dir_out1, aliquot_tmp, ".lymphoid_reclustered.", "dotplot.",run_id, ".png")
      png(file = file2write, width = 4000, height = 1200, res = 150)
      print(p)
      dev.off()
    }
  }
}

