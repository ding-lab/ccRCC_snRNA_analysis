# Yige Wu @WashU Oct 2019
## for integrating two snRNA datasets for sample CPT0086820004 and CPT0075130004 (from cellranger output with premrna reference)

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
gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Kidney_Markers/RCC_marker_gene_table_and_literature_review - Gene2CellType_Tab_w.HumanProteinAtlas.20191022.v4.tsv")

# input non-immune integrated object --------------------------------------
object2plot <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/process_nonimmune_cells_on_cluster/20191022.v1/NonImmune_Integrated.20191022.v1.RDS")

# input DEG ---------------------------------------------------------------
deg_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/process_nonimmune_cells_on_cluster/20191022.v1/NonImmune.DEGs.Pos.txt", data.table = F)
deg_tab$gene %>% head()

# input pathway annotated table -------------------------------------------
reactomepa_sup_table <- NULL
for (cluster_tmp in cluster2celltype_tab$Cluster) {
  reactomepa_tab_tmp <- fread(input = paste0(dir_out, "Cluster", cluster_tmp, "_", cluster2celltype_tab$Enriched_Cell_Type_Abbr[cluster2celltype_tab$Cluster == cluster_tmp], "_DEG_Enriched_Pathways.tsv"))
  reactomepa_tab_tmp$cluster <-cluster_tmp 
  reactomepa_sup_table <- rbind(reactomepa_sup_table, reactomepa_tab_tmp)
}

# re-annotate the genes to have only one pathway annotated to one gene --------
gene_tmp <- "SLC6A13"
reactomepa_gene_tmp <- reactomepa_sup_table %>%
  filter(grepl(pattern = gene_tmp, x = geneID)) %>%
  arrange(desc(Count))

deg_tab$reactome_pathway <- sapply(deg_tab$gene, function(g, pathway_tab) {
  reactomepa_gene_tmp <- pathway_tab %>%
    filter(grepl(pattern = g, x = geneID)) %>%
    arrange(desc(Count))
  if (nrow(reactomepa_gene_tmp) == 0) {
    return("Other")
  } else {
    return(reactomepa_gene_tmp$Description[1])
  }
}, pathway_tab = reactomepa_sup_table)


# filter DEG --------------------------------------------------------------



# dotplot for tumor only -----------------------------------------------------------------
DefaultAssay(object2plot) <- "RNA"
object2plot_tmp <- subset(x = object2plot, subset = seurat_clusters %in% cluster2celltype_tab$Cluster[grepl(x = cluster2celltype_tab$Enriched_Cell_Type_Abbr, pattern = "ccRCC")])
gene2p_tab <- deg_tab %>%
  filter(p_val_adj < 0.05) %>%
  filter(reactome_pathway != "Other")

gene2p_tab <- gene2p_tab %>%
  group_by(cluster) %>%
  dplyr::top_n(n = 10, wt = avg_logFC)

genes2plot <- unique(gene2p_tab$gene)

p <- DotPlot(object = object2plot_tmp, features = genes2plot, col.min = 0)
p$data$gene_pathway<- plyr::mapvalues(p$data$features.plot, from = deg_tab$gene, to = deg_tab$reactome_pathway)
p$data$cluster_cell_type <- plyr::mapvalues(p$data$id, from = cluster2celltype_tab$Cluster, to = cluster2celltype_tab$Enriched_Cell_Type_Abbr)
p <- p  + RotatedAxis()
p <- p + facet_grid(cluster_cell_type~gene_pathway, scales = "free", space = "free", drop = T)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               strip.text.x = element_text(angle = 90, vjust = 0.5, face = "bold"),
               strip.text.y = element_text(angle = 0, vjust = 0.5),
               axis.text.x = element_text(size = 14, face = "bold"),
               strip.placement = "outside")
file2write <- paste0(dir_out, "Dotplot_Tumor_Top10DEG_Exp.", run_id, ".png")
png(file = file2write, width = 5000, height = 2000, res = 150)
print(p)
dev.off()


# dotplot for nonimmune -----------------------------------------------------------------
DefaultAssay(object2plot) <- "RNA"

gene2p_tab <- deg_tab %>%
  filter(p_val_adj < 0.05) %>%
  filter(reactome_pathway != "Other")

gene2p_tab <- gene2p_tab %>%
  group_by(cluster) %>%
  dplyr::top_n(n = 15, wt = avg_logFC)

genes2plot <- unique(gene2p_tab$gene)

p <- DotPlot(object = object2plot, features = genes2plot, col.min = 0)
p$data$gene_pathway<- plyr::mapvalues(p$data$features.plot, from = deg_tab$gene, to = deg_tab$reactome_pathway)
p$data$cluster_cell_type <- plyr::mapvalues(p$data$id, from = cluster2celltype_tab$Cluster, to = cluster2celltype_tab$Enriched_Cell_Type_Abbr)
p <- p  + RotatedAxis()
p <- p + facet_grid(cluster_cell_type~gene_pathway, scales = "free", space = "free", drop = T)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               strip.text.x = element_text(angle = 90, vjust = 0.5),
               strip.text.y = element_text(angle = 0, vjust = 0.5),
               axis.text.x = element_text(size = 9),
               strip.placement = "outside")
file2write <- paste0(dir_out, "Dotplot_NonImmune_Marker_Exp.", run_id, ".png")
png(file = file2write, width = 5000, height = 2000, res = 150)
print(p)
dev.off()


