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


# input genes to plot -----------------------------------------------------
IAP_genes <- c("BIRC2", "BIRC3", "XIAP", "BIRC7", "BIRC5")
cdk_genes <- c("CDK1", "CDK2", "CDK4", "CDK5", "CDK6", "CDK7", "CDK9", "WEE1")
mtor_genes <- c("MTOR", "AKTS1", "DEPTOR", "RPTOR", "MLST8", "MAPKAP1", "PRR5", "RICTOR")
## targets by carbozantinib: https://www.drugbank.ca/drugs/DB08875
rtk_genes <- c("MET", "AXL", "FLT1", "KDR", "FLT3", "KIT", "RET", "NTRK2", "TEK")
## specify the genes to be plotted 
drug_targets <- c(IAP_genes, cdk_genes, mtor_genes, rtk_genes)


# generate a table with average expression for all HIF targets for each cluster------------
p <- DotPlot(object = object2plot, features = drug_targets, col.min = 0)
exp_stat_tab <- p$data
# write.table(x = exp_stat_tab, paste0(dir_out, "HIF_Targets_snRNA_Exp_By_Cluster.Integration." ,integration_id, ".Run.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)

# plot dotplot for only differentially expressed HIF targets --------------
genes2plot <- intersect(drug_targets, unique(deg_tab$gene))
genes2plot <- drug_targets

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

file2write <- paste0(dir_out, "Dotplot_Emerson_Drug_Targets_snRNA.", run_id, ".pdf")
pdf(file = file2write, width = 14, height = 8, useDingbats = F)
print(p)
dev.off()

file2write <- paste0(dir_out, "Dotplot_Emerson_Drug_Targets_snRNA.", run_id, ".png")
png(file = file2write, width = 2100, height = 1200, res = 150)
print(p)
dev.off()


