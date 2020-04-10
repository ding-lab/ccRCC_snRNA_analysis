# Yige Wu @WashU Nov 2019
## for plotting the marker genes for individual sample in dotplot

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
# snRNA_aliquot_ids <- c("CPT0019130004", "CPT0001260013", "CPT0086350004", "CPT0010110013", "CPT0001180011", "CPT0025890002", "CPT0075140002", "CPT0020120013", "CPT0001220012", "CPT0014450005")
# snRNA_aliquot_ids <- c("CPT0025880013", "CPT0015810004", "CPT0075130004", "CPT0086820004")
snRNA_aliquot_ids <- c("CPT0010110013")


# input cluster 2 cell type tables ----------------------------------------
cluster2celltype_allcells_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - Individual.AllCluster2Cell_Type.20191126.v1.tsv", data.table = F)

# input ligand-receptor pairs ---------------------------------------------
lig_rep_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Ligand_Receptor/PairsLigRec.txt", data.table = F)


# set recluster_tumor_id --------------------------------------------------
recluster_tumor_id <- "20191119.v1"

# input seurat processing summary ------------------------------------------------
seurat_summary <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/summary/ccRCC_snRNA_Downstream_Processing - Seurat_Preprocessing.20191125.v1.tsv", data.table = F)
seurat_summary2process <- seurat_summary %>%
  filter(Cellranger_reference_type == "pre-mRNA") %>%
  filter(Proceed_for_downstream == "Yes") %>%
  filter(Aliquot %in% snRNA_aliquot_ids) %>%
  mutate(Path_seurat_object = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                     "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                     "/", Aliquot, FACS, "_processed.rds")) %>%
  mutate(Path_seurat_object_malignant = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/individual_cluster/recluster_tumor/", recluster_tumor_id, "/", 
                                               Aliquot, FACS, "/", Aliquot, FACS, ".Malignant_Reclustered.", recluster_tumor_id, ".RDS")) %>%
  mutate(Path_deg_table = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                 "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                 "/", Aliquot, FACS, ".DEGs.Pos.txt"))
seurat_summary2process$Path_seurat_object

# plot dotplot by sample --------------------------------------------------
for (snRNA_aliquot_id_tmp in snRNA_aliquot_ids) {
  ## input seurat object
  seurat_obj_path <- seurat_summary2process$Path_seurat_object[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
  seurat_obj <- readRDS(file = seurat_obj_path)
  DefaultAssay(seurat_obj) <- "RNA"
  
  ## input DEG
  deg_tab_path <- seurat_summary2process$Path_deg_table[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
  deg_tab_path
  deg_tab <- fread(input = deg_tab_path, data.table = F)
  deg_tab <- deg_tab %>%
    filter(p_val_adj < 0.05) %>%
    filter(avg_logFC > 0.5)
  
  ## get differentially expression genes
  deg_genes <- unique(deg_tab$gene)
  
  ## plot marker gene expression
  lig_rep_tab2plot <- lig_rep_tab %>%
    filter(Ligand.ApprovedSymbol %in% deg_genes) %>%
    filter(Receptor.ApprovedSymbol %in% deg_genes)
  
  ## divide cell group
  ### group stromal and immune cells
  cluster2celltype_allcells_tab_aliquot <- cluster2celltype_allcells_tab %>%
    filter(Aliquot == snRNA_aliquot_id_tmp)
  
  ### get the cluster numbers by cell type category
  malignant_clusters <- cluster2celltype_allcells_tab_aliquot$Cluster[cluster2celltype_allcells_tab_aliquot$Is_Malignant == "Yes"]
  stromal_clusters <- cluster2celltype_allcells_tab_aliquot$Cluster[cluster2celltype_allcells_tab_aliquot$Is_Stromal == "Yes"]
  immune_clusters <- cluster2celltype_allcells_tab_aliquot$Cluster[cluster2celltype_allcells_tab_aliquot$Is_Immune == "Yes"]
  
  ### create text for cell type
  cell2celltype_tab <- data.frame(barcode = rownames(seurat_obj@meta.data), allcells_cluster = seurat_obj@meta.data$seurat_clusters)
  cell2celltype_tab$celltype_text <- mapvalues(x = cell2celltype_tab$allcells_cluster, from = cluster2celltype_allcells_tab_aliquot$Cluster, to = cluster2celltype_allcells_tab_aliquot$Cell_Type_Abbr)
  cell2celltype_tab$celltype_text <- paste0(cell2celltype_tab$celltype_text, "_C", cell2celltype_tab$allcells_cluster)
  
  
  ### group malignant cells into new tumor clusters
  #### input re-intergrated object for tumor cells
  malignant_seurat_obj_path <- seurat_summary2process$Path_seurat_object_malignant[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
  malignant_seurat_obj <- readRDS(file = malignant_seurat_obj_path)
  
  #### group maligant cells into new tumor groups
  malignant_barcodes <- cell2celltype_tab$barcode[cell2celltype_tab$allcells_cluster %in% malignant_clusters]
  malignant_new_clusters <- mapvalues(x = malignant_barcodes, from = rownames(malignant_seurat_obj@meta.data), to = as.vector(malignant_seurat_obj@meta.data$seurat_clusters))
  cell2celltype_tab$celltype_text[cell2celltype_tab$barcode %in% malignant_barcodes] <- paste0("MalignantPT", "_C", malignant_new_clusters)
  
  ## add celltype category to the cell2celltype table
  cell2celltype_tab$celltype_cat <- ifelse(cell2celltype_tab$allcells_cluster %in% malignant_clusters, "Malignant",
                                           ifelse(cell2celltype_tab$allcells_cluster %in% stromal_clusters, "Stromal",
                                                  ifelse(cell2celltype_tab$allcells_cluster %in% immune_clusters, "Immune", "Other")))
  
  ## add celltype label to the seurat object
  seurat_obj@meta.data$celltype_text <- mapvalues(x = rownames(seurat_obj@meta.data), from = cell2celltype_tab$barcode, to = cell2celltype_tab$celltype_text)
  
  ## get the DEG genes by cell type categories
  deg_tab_uniq <- deg_tab %>%
    arrange(gene, -avg_logFC)
  deg_tab_uniq <- deg_tab_uniq[!duplicated(deg_tab_uniq$gene),]
  
  malignant_deg_genes <- unique(deg_tab_uniq$gene[deg_tab_uniq$cluster %in% malignant_clusters])
  stromal_deg_genes <- unique(deg_tab_uniq$gene[deg_tab_uniq$cluster %in% stromal_clusters])
  immune_deg_genes <- unique(deg_tab_uniq$gene[deg_tab_uniq$cluster %in% immune_clusters])
  
  ## draw expression for the ligands
  ### define the ligand genes to plot
  ligands2plot <- unique(lig_rep_tab2plot$Ligand.ApprovedSymbol)
  ### define the cell type the genes are most highly expressed
  ligand_celltype_exp <- ifelse(ligands2plot %in% malignant_deg_genes, "Malignant_Expressed",
                                ifelse(ligands2plot %in% stromal_deg_genes, "Stromal_Expressed", "Immune_Expressed"))
  ### make plots
  p <- DotPlot(object = seurat_obj, features = ligands2plot, col.min = 0, group.by = "celltype_text")
  p$data$celltype_cat <- mapvalues(x = p$data$id, from = cell2celltype_tab$celltype_text, to = cell2celltype_tab$celltype_cat)
  p$data$gene_celltype_exp <- mapvalues(x = p$data$features.plot, from = ligands2plot, to = ligand_celltype_exp)
  p$data$gene_celltype_exp <- factor(p$data$gene_celltype_exp, levels = c("Malignant_Expressed", "Stromal_Expressed", "Immune_Expressed"))
  p <- p + facet_grid(celltype_cat~gene_celltype_exp, scales = "free", space = "free", drop = T, switch = "y")
  p <- p + theme(panel.spacing = unit(0, "lines"),
                 strip.background.y = element_rect(colour = "black", fill = "white"),
                 strip.background.x = element_rect(colour = "black", fill = "white"),
                 panel.border = element_rect(colour = "black"),
                 strip.text.y = element_text(angle = 180, vjust = 0.5),
                 strip.placement = "outside")
  p <- p + theme(axis.text.x = element_text(angle = 90, size = 12),
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank())
  p <- p + theme(panel.grid.major = element_line(colour = "black"))
  p <- p + theme(legend.position = "top")
  ligand_exp_plot <- p
  
  ## draw expression for the receptors
  ### define the receptor genes to plot
  receptors2plot <- unique(lig_rep_tab2plot$Receptor.ApprovedSymbol)
  ### define the cell type the genes are most highly expressed
  receptor_celltype_exp <- ifelse(receptors2plot %in% malignant_deg_genes, "Malignant_Expressed",
                                  ifelse(receptors2plot %in% stromal_deg_genes, "Stromal_Expressed", "Immune_Expressed"))
  ### make plots
  p <- DotPlot(object = seurat_obj, features = receptors2plot, col.min = 0, group.by = "celltype_text")
  p$data$celltype_cat <- mapvalues(x = p$data$id, from = cell2celltype_tab$celltype_text, to = cell2celltype_tab$celltype_cat)
  p$data$gene_celltype_exp <- mapvalues(x = p$data$features.plot, from = receptors2plot, to = receptor_celltype_exp)
  p$data$gene_celltype_exp <- factor(p$data$gene_celltype_exp, levels = c("Malignant_Expressed", "Stromal_Expressed", "Immune_Expressed"))
  p <- p + facet_grid(celltype_cat~gene_celltype_exp, scales = "free", space = "free", drop = T, switch = "y")
  p <- p + theme(panel.spacing = unit(0, "lines"),
                 strip.background.y = element_rect(colour = "black", fill = "white"),
                 strip.background.x = element_rect(colour = "black", fill = "white"),
                 panel.border = element_rect(colour = "black"),
                 strip.text.y = element_text(angle = 180, vjust = 0.5),
                 strip.placement = "outside")
  p <- p + theme(axis.text.x = element_text(angle = 90, size = 12),
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank())
  p <- p + theme(panel.grid.major = element_line(colour = "black"))
  p <- p + theme(legend.position = "bottom")
  receptor_exp_plot <- p
  receptor_exp_plot
  
  ## draw lines between ligand and receptors
  ### get the ordered list of ligand genes in the ligand plot
  ligands2plot_ordered <- c(ggplot_build(ligand_exp_plot)$layout$panel_params[[1]]$x.labels,
                            ggplot_build(ligand_exp_plot)$layout$panel_params[[2]]$x.labels,
                            ggplot_build(ligand_exp_plot)$layout$panel_params[[3]]$x.labels)
  
  receptors2plot_ordered <- c(ggplot_build(receptor_exp_plot)$layout$panel_params[[1]]$x.labels,
                              ggplot_build(receptor_exp_plot)$layout$panel_params[[2]]$x.labels,
                              ggplot_build(receptor_exp_plot)$layout$panel_params[[3]]$x.labels)
  
  
  ### annotate the pair type by who is expressing the ligand and who is expressing the receptor
  lig_rep_tab2plot$pair_by_celltype_exp <- paste0(mapvalues(x = lig_rep_tab2plot$Ligand.ApprovedSymbol, from = ligands2plot, to = gsub(x = ligand_celltype_exp, pattern = "_Expressed", replacement = "")),
                                                  "->",
                                                  mapvalues(x = lig_rep_tab2plot$Receptor.ApprovedSymbol, from = receptors2plot, to = gsub(x = receptor_celltype_exp, pattern = "_Expressed", replacement = "")))
  
  ### set color responding to each type
  line_colors <- brewer.pal(9,"Set1")
  names(line_colors) <- unique(lig_rep_tab2plot$pair_by_celltype_exp)
  
  line_plot_tab <- lig_rep_tab2plot %>%
    select(Ligand.ApprovedSymbol, Receptor.ApprovedSymbol, Pair.Name, pair_by_celltype_exp)
  line_plot_tab <- melt(line_plot_tab, id.vars = c("Pair.Name", "pair_by_celltype_exp"))
  line_plot_tab$x <- sapply(X = line_plot_tab$value, FUN = function(g, ligand_list, receptor_list) {
    if (g %in% ligands2plot_ordered) {
      x <- which(ligand_list == g)
      return(x)
    } else {
      x <- (which(receptor_list == g)/length(receptor_list))*length(ligand_list)
      return(x)
    }}, ligand_list = ligands2plot_ordered, receptor_list = receptors2plot_ordered)
  line_plot_tab <- line_plot_tab %>%
    mutate(gene_type = str_split_fixed(string = variable, pattern = ".ApprovedSymbol", n = 2)[,1])
  line_plot_tab$y <- ifelse(line_plot_tab$gene_type == "Ligand", 1, 0) 
  
  p <- ggplot()
  p <- p + geom_line(data = line_plot_tab, mapping = aes(x = x, y = y, group = Pair.Name, color = pair_by_celltype_exp))
  p <- p + scale_color_manual(values = line_colors)
  p <- p + geom_text(data = line_plot_tab[line_plot_tab$gene_type == "Ligand",], mapping = aes(x = x, y = y, label = value), angle = 90, nudge_y = 0.05)
  p <- p + geom_text(data = line_plot_tab[line_plot_tab$gene_type == "Receptor",], mapping = aes(x = x, y = y, label = value), angle = 90, nudge_y = -0.05)
  p <- p + scale_x_continuous(expand = c(0.01,0.01), limits=range(line_plot_tab$x))
  p <- p + theme_bw()
  p <- p + theme(axis.title = element_blank(),
                 axis.ticks = element_blank(),
                 axis.text = element_blank(),
                 panel.grid.major = element_blank())
  p <- p + theme(legend.position = c(-0.1,0.5),
                 legend.text = element_text(size = 15, face = "bold"))
  line_plot <- p
  
  file2write <- paste(dir_out, snRNA_aliquot_id_tmp,".Individual_Clustered.", "LigandReceptor.Dotplot.", run_id, ".png", sep="")
  png(file2write, width = 2500, height = 3000, res = 150)
  p <- plot_grid(ligand_exp_plot, line_plot, receptor_exp_plot, 
                 nrow = 3, rel_heights = c(1, 1, 1), align = "v", axis = "lr", greedy = F)
  print(p)
  dev.off()
  
  for (pair_type_tmp in unique(lig_rep_tab2plot$pair_by_celltype_exp)) {
    p <- ggplot()
    p <- p + geom_line(data = line_plot_tab[line_plot_tab$pair_by_celltype_exp == pair_type_tmp,], mapping = aes(x = x, y = y, group = Pair.Name), color = line_colors[pair_type_tmp])
    p <- p + geom_text(data = line_plot_tab[line_plot_tab$gene_type == "Ligand",], mapping = aes(x = x, y = y, label = value), angle = 90, nudge_y = 0.05)
    p <- p + geom_text(data = line_plot_tab[line_plot_tab$gene_type == "Receptor",], mapping = aes(x = x, y = y, label = value), angle = 90, nudge_y = -0.05)
    p <- p + scale_x_continuous(expand = c(0.01,0.01), limits=range(line_plot_tab$x))
    p <- p + theme_bw()
    p <- p + theme(axis.title = element_blank(),
                   axis.ticks = element_blank(),
                   axis.text = element_blank(),
                   panel.grid.major = element_blank())
    p <- p + theme(legend.position = c(-0.1,0.5),
                   legend.text = element_text(size = 15, face = "bold"))
    line_plot_tmp <- p
    
    file2write <- paste(dir_out, snRNA_aliquot_id_tmp,".Individual_Clustered.", "LigandReceptor.", pair_type_tmp, ".Dotplot.", run_id, ".png", sep="")
    png(file2write, width = 2500, height = 3000, res = 150)
    p <- plot_grid(ligand_exp_plot, line_plot_tmp, receptor_exp_plot, 
                   nrow = 3, rel_heights = c(1, 1, 1), align = "v", axis = "lr", greedy = F)
    print(p)
    dev.off()
  }
}


