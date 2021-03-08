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

# input non-immune integrated object --------------------------------------
object2plot <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/integrate_seurat_objects/20191021.v1/Renal_Integrated.20191021.v1.RDS")


# input cluster 2 cell type table -----------------------------------------
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - AllCluster2Cell_Type.20191022.v1.tsv", data.table = F)


# annotate cell type for each barcode -------------------------------------
object2plot@meta.data$cluster_cell_type <- plyr::mapvalues(object2plot@meta.data$seurat_clusters, from = cluster2celltype_tab$Cluster, to = cluster2celltype_tab$Enriched_Cell_Type_Abbr)

# Dimplot split by aliquot ------------------------------------------------
p <- DimPlot(object2plot, reduction = "umap", group.by = "cluster_cell_type", label = T)
file2write <- paste(dir_out, "Dimplot_by_Clusters_with_CellType.", run_id, ".png", sep="")
png(file = file2write, width = 2000, height = 1000, res = 150)
print(p)
dev.off()

# Dimplot split by aliquot ------------------------------------------------
p <- DimPlot(object2plot, reduction = "umap", group.by = "cluster_cell_type", label = T, label.size	= 5, repel = T, 
             split.by = "orig.ident", ncol = 4)
file2write <- paste(dir_out, "Dimplot_by_Clusters_and_Aliquots_with_CellType.4col.", run_id, ".png", sep="")
png(file = file2write, width = 6000, height = 3000, res = 150)
print(p)
dev.off()
