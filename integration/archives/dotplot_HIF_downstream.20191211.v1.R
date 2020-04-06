# Yige Wu @WashU Oct 2019
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
gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Kidney_Markers/RCC_marker_gene_table_and_literature_review - Gene2CellType_Tab_w.HumanProteinAtlas.20191022.v5.tsv")

# input non-immune integrated object --------------------------------------
integration_id <- "20191021.v1"
object2plot <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/integrate_seurat_objects/20191021.v1/Renal_Integrated.20191021.v1.RDS")

# input DEG ---------------------------------------------------------------
deg_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/integrate_seurat_objects/20191021.v1/Renal.DEGs.Pos.txt", data.table = F)

# input cluster 2 cell type table -----------------------------------------
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - Intergrated.AllCluster2Cell_Type.20191211.v1.tsv", data.table = F)

# specify the HIF downstream to plot ----------------------------------------
tf_tab <- fread(input = "./Ding_Lab/Projects_Current/TP53_shared_data/resources/PPI/TF_interactions.txt", data.table = F)
HIF_tf_tab <- tf_tab %>%
  filter(source_genesymbol %in% c("HIF1A", "EPAS1"))
HIF_targets <- HIF_tf_tab$target_genesymbol
HIF_targets <- unique(c("PLIN2", "CA12", "PLG", "IL6", "ADM", "VEGFA", "BNIP3", "HK1", "HK2", "PFK", "ALDOA", "PGK1", "LDHA", "NOS2", "ABL2", "EPO", "POUSF1", "SCGB3A1", "TGFA", "CCND1", "DLL4", "ANGPT2",
                        HIF_targets))

# generate a table with average expression for all HIF targets for each cluster------------
p <- DotPlot(object = object2plot, features = HIF_targets, col.min = 0)
exp_stat_tab <- p$data
write.table(x = exp_stat_tab, paste0(dir_out, "HIF_Targets_snRNA_Exp_By_Cluster.Integration." ,integration_id, ".Run.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)

# plot dotplot for only differentially expressed HIF targets --------------
genes2plot <- intersect(HIF_targets, unique(deg_tab$gene))

## divide cell group
### get the cluster numbers by cell type category
malignant_clusters <- cluster2celltype_tab$Cluster[cluster2celltype_tab$Is_Malignant == "Yes"]
stromal_clusters <- cluster2celltype_tab$Cluster[cluster2celltype_tab$Is_Stromal == "Yes"]
immune_clusters <- cluster2celltype_tab$Cluster[cluster2celltype_tab$Is_Immune == "Yes"]

## add celltype category to the cell2celltype table
cluster2celltype_tab$celltype_cat <- ifelse(cluster2celltype_tab$Cluster %in% malignant_clusters, "Malignant",
                                         ifelse(cluster2celltype_tab$Cluster %in% stromal_clusters, "Stromal",
                                                ifelse(cluster2celltype_tab$Cluster %in% immune_clusters, "Immune", "Other")))

## get the DEG genes by cell type categories
deg_tab_uniq <- deg_tab %>%
  arrange(gene, -avg_logFC)
deg_tab_uniq <- deg_tab_uniq[!duplicated(deg_tab_uniq$gene),]

malignant_deg_genes <- unique(deg_tab_uniq$gene[deg_tab_uniq$cluster %in% malignant_clusters])
stromal_deg_genes <- unique(deg_tab_uniq$gene[deg_tab_uniq$cluster %in% stromal_clusters])
immune_deg_genes <- unique(deg_tab_uniq$gene[deg_tab_uniq$cluster %in% immune_clusters])

### define the cell type the genes are most highly expressed
gene_celltype_exp_cats <- ifelse(genes2plot %in% malignant_deg_genes, "Malignant_Expressed",
                              ifelse(genes2plot %in% stromal_deg_genes, "Stromal_Expressed", "Immune_Expressed"))

## actual plotting
DefaultAssay(object2plot) <- "RNA"
p <- DotPlot(object = object2plot, features = genes2plot, col.min = 0)
# p$data$HIF1A_downstream <- ifelse(p$data$features.plot %in% HIF_tf_tab$target_genesymbol[HIF_tf_tab$source_genesymbol == "HIF1A" & HIF_tf_tab$is_stimulation == 1], "HIF1A_Stimulated", 
#                                   ifelse(p$data$features.plot %in% HIF_tf_tab$target_genesymbol[HIF_tf_tab$source_genesymbol == "HIF1A" & HIF_tf_tab$is_inhibition == 1], "HIF1A_Inhibited", "HIF1A_Regulated"))
# p$data$EPAS1_downstream <- ifelse(p$data$features.plot %in% HIF_tf_tab$target_genesymbol[HIF_tf_tab$source_genesymbol == "EPAS1" & HIF_tf_tab$is_stimulation == 1], "EPAS1_Stimulated", 
#                                   ifelse(p$data$features.plot %in% HIF_tf_tab$target_genesymbol[HIF_tf_tab$source_genesymbol == "EPAS1" & HIF_tf_tab$is_inhibition == 1], "EPAS1_Inhibited", "EPAS1_Regulated"))

p$data$celltype_text <- plyr::mapvalues(p$data$id, from = cluster2celltype_tab$Cluster, to = cluster2celltype_tab$Enriched_Cell_Type_Abbr)
p$data$celltype_cat <- mapvalues(x = p$data$id, from = cluster2celltype_tab$Cluster, to = cluster2celltype_tab$celltype_cat)
p$data$gene_celltype_exp_cat <- mapvalues(x = p$data$features.plot, from = genes2plot, to = gene_celltype_exp_cats)
p$data$gene_celltype_exp_cat <- factor(p$data$gene_celltype_exp_cat, levels = c("Malignant_Expressed", "Stromal_Expressed", "Immune_Expressed"))
p <- p  + RotatedAxis()
p <- p + facet_grid(celltype_cat + celltype_text~gene_celltype_exp_cat, scales = "free", space = "free", drop = T)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background.y = element_rect(colour = "black", fill = "white"),
               strip.background.x = element_rect(colour = "black", fill = "white"),
               panel.border = element_rect(colour = "black"),
               strip.text.x = element_text(angle = 90, vjust = 0.5),
               strip.text.y = element_text(angle = 0, vjust = 0.5),
               axis.text.x = element_text(size = 12, face = "bold"),
               strip.placement = "outside")
file2write <- paste0(dir_out, "Dotplot_HIF_Downstream_Exp.", run_id, ".png")
png(file = file2write, width = 2500, height = 1300, res = 150)
print(p)
dev.off()


