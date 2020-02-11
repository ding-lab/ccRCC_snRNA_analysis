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
object2plot <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/integrate_seurat_objects/20191021.v1/Renal_Integrated.20191021.v1.RDS")

# input DEG ---------------------------------------------------------------
deg_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/integrate_seurat_objects/20191021.v1/Renal.DEGs.Pos.txt", data.table = F)

# input cluster 2 cell type table -----------------------------------------
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - AllCluster2Cell_Type.20191022.v1.tsv", data.table = F)
nonimmune_clusters <- cluster2celltype_tab %>%
  filter(Is_Immune == "No") %>%
  select(Cluster)
nonimmune_clusters <- nonimmune_clusters$Cluster


# input TF table ----------------------------------------------------------
tf_tab <- fread(input = "./Ding_Lab/Projects_Current/TP53_shared_data/resources/PPI/TF_interactions.txt", data.table = F)
HIF_tf_tab <- tf_tab %>%
  filter(source_genesymbol %in% c("HIF1A", "EPAS1"))

# Dotplot for all druggable targets -------------------------------------------
DefaultAssay(object2plot) <- "RNA"

genes2plot <- unique(HIF_tf_tab$target_genesymbol)

p <- DotPlot(object = object2plot, features = genes2plot, col.min = 0)
p$data$HIF1A_downstream <- ifelse(p$data$features.plot %in% HIF_tf_tab$target_genesymbol[HIF_tf_tab$source_genesymbol == "HIF1A" & HIF_tf_tab$is_stimulation == 1], "HIF1A_Stimulated", 
                                  ifelse(p$data$features.plot %in% HIF_tf_tab$target_genesymbol[HIF_tf_tab$source_genesymbol == "HIF1A" & HIF_tf_tab$is_inhibition == 1], "HIF1A_Inhibited", "HIF1A_Regulated"))
p$data$EPAS1_downstream <- ifelse(p$data$features.plot %in% HIF_tf_tab$target_genesymbol[HIF_tf_tab$source_genesymbol == "EPAS1" & HIF_tf_tab$is_stimulation == 1], "EPAS1_Stimulated", 
                                  ifelse(p$data$features.plot %in% HIF_tf_tab$target_genesymbol[HIF_tf_tab$source_genesymbol == "EPAS1" & HIF_tf_tab$is_inhibition == 1], "EPAS1_Inhibited", "EPAS1_Regulated"))

# p$data$HIF1A_downstream <- plyr::mapvalues(p$data$features.plot, from = HIF_tf_tab$target_genesymbol, to = gene2cellType_tab$Biomarker_Type)
p$data$cluster_cell_type <- plyr::mapvalues(p$data$id, from = cluster2celltype_tab$Cluster, to = cluster2celltype_tab$Enriched_Cell_Type_Abbr)
p$data$is_immune <- ifelse(p$data$id %in% nonimmune_clusters, "Non_Immune", "Immune")
p$data$is_immune <- factor(p$data$is_immune, levels = c("Non_Immune", "Immune"))

p <- p  + RotatedAxis()
p <- p + facet_grid(is_immune + cluster_cell_type~HIF1A_downstream + EPAS1_downstream, scales = "free", space = "free", drop = T)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               strip.text.x = element_text(angle = 90, vjust = 0.5),
               strip.text.y = element_text(angle = 0, vjust = 0.5),
               axis.text.x = element_text(size = 9),
               strip.placement = "outside")

file2write <- paste0(dir_out, "Dotplot_HIF_Downstream_Exp.", run_id, ".png")
png(file = file2write, width = 2500, height = 1300, res = 150)
print(p)
dev.off()


