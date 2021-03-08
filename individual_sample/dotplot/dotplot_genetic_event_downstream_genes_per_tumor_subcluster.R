# Yige Wu @WashU Feb 2020
## for plotting the expression of the potential downstream genes in tumor cells

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
## input seurat object paths
srat_paths <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/recluster_nephron_epithelium_cells_in_individual_samples/20200225.v1/Seurat_Object_Paths.Malignant_Nephron_Epithelium20200225.v1.tsv", data.table = F)
## input barcode to cluster info
barcode2celltype_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/fetch_data_malignant_nephron_epithelium_cells_reclustered/20200225.v1/Cluster_UMAP_Data.Malignant_Nephron_Epithelium.20200225.v1.tsv", data.table = F)
## input the gene list to plot
genetic_alt_downstream_genes <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/dependencies/write_ccrcc_genetic_event_downstream_genes/20200227.v1/ccRCC_Genetic_Event_Downstream_Genes.20200227.v1.tsv", data.table = F)
## input per case bulk CNV profile
bulk_cnv_state_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/bulk/write_sample_cnv_profile/20200227.v1/Bulk_WGS_Chr_CNV_Profile.20200227.v1.tsv", data.table = F)
## set the mininal pct expressed in any cluster
min.exp.pct <- 25

# get expression data for the gene list -----------------------------------
aliquot_tmp <- "CPT0000880001"
for (aliquot_tmp in srat_paths$Aliquot) {
  ## input individually processed seurat object
  seurat_obj_path <- srat_paths$Path_seurat_object[srat_paths$Aliquot == aliquot_tmp]
  seurat_obj_path
  seurat_object <- readRDS(file = seurat_obj_path)
  
  ## input differentially expressed genes across tumor subclusters
  deg_df <- fread(input = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/findallmarkers_malignant_nephron_epithelium_cells_reclustered/20200225.v1/", aliquot_tmp, ".FindAllMarkers.Wilcox.Pos.tsv"))
  deg_df <- deg_df %>%
    filter(p_val_adj < 0.05)
  
  ## filter genes to those within DEGs
  genes2plot_filtered <- unique(genetic_alt_downstream_genes$target_genesymbol)
  genes2plot_filtered <- genes2plot_filtered[genes2plot_filtered %in% deg_df$gene]
  
  ## plot who is whose target
  plot_data_df2 <- genetic_alt_downstream_genes
  ### only show the 
  plot_data_df2 <- plot_data_df2 %>%
    filter(target_genesymbol %in% genes2plot_filtered)
  ### reformat the long data frame to wide data frame to get the TFs for each target gene
  plot_data_mat2 <- dcast(data = plot_data_df2, formula = source_genesymbol ~ target_genesymbol, fill = "", value.var = "source_genesymbol")
  plot_data_mat2 <- plot_data_mat2[,-1]
  target2sources <- sapply(colnames(plot_data_mat2), function(col_name, wide_df) {
    text_tmp <- wide_df[, col_name]
    text_tmp <- text_tmp[text_tmp != ""]
    text_tmp <- paste0(text_tmp, collapse = "|")
    return(text_tmp)
  }, wide_df = plot_data_mat2)
  target2sources_df <- data.frame(target_genesymbol = names(target2sources),
                                  source_genesymbols = target2sources)
  target2sources_df <- target2sources_df %>%
    arrange(source_genesymbols)
  ### order the target gene symbols
  plot_data_df2$target_genesymbol <- factor(plot_data_df2$target_genesymbol, levels = target2sources_df$target_genesymbol)
  ### ggplot
  p <- ggplot()
  p <- p + geom_tile(data = plot_data_df2, mapping = aes(x = target_genesymbol, y = source_genesymbol))
  p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 12, face = "bold"),
                 axis.text.y = element_text(size = 12, face = "bold"))
  p2 <- p
  p2
  
  ## dotplot
  p <- DotPlot(object = seurat_object, features = genes2plot_filtered)
  p$data$source_genesymbols <- mapvalues(x = p$data$features.plot, from = target2sources_df$target_genesymbol, to = as.vector(target2sources_df$source_genesymbols))
  p$data$features.plot <- factor(p$data$features.plot, levels = target2sources_df$target_genesymbol)
  p <- p + ggtitle(label = paste0("", aliquot_tmp, " Tumor Subclusters"), subtitle = paste0("CNV Downstream Gene Expression"))
  p <- p + facet_grid(.~source_genesymbols, scales = "free", space = "free", drop = T)
  p <- p + theme_bw()
  p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 12, face = "bold"),
                 axis.text.y = element_text(size = 12, face = "bold"))
  p <- p + theme(panel.spacing = unit(0, "lines"))
  p <- p + theme(strip.text.y = element_text(angle = 0))
  p
  p1 <- p
  p1
  
  stop("")
}

