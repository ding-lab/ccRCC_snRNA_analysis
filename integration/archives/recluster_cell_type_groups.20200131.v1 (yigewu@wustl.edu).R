# Yige Wu @WashU Jan 2020
## Re-cluster each cell type group for the tumor-normal integrated dataset

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

# set case id -------------------------------------------------------------
case_id_tmp <- "C3N-01200"

# input seurat object -----------------------------------------------------
seurat_obj_path <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/integrate_same_patient_tumor_normal/20200117.v1/", case_id_tmp, ".Tummor_Normal.Integrated.20200117.v1.RDS")
## input seurat object
orig.object <- readRDS(file = seurat_obj_path)
DefaultAssay(orig.object) <- "integrated"

# get clusters to be processed --------------------------------------------
## input cluster to cell type table
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - C3N-01200.TumorNormal.Integrated.20200131.v1.tsv", data.table = F)

## set cell type group name to be processed
cell_type_group_tmp <- "Nephron_Epithelium"

## get clusters to be processed by cell type group
clusters2process <- cluster2celltype_tab[cluster2celltype_tab[,cell_type_group_tmp] == "Yes", "Cluster"] 
clusters2process

# subset object by non-immune clusters ------------------------------------
new.object <- subset(orig.object, idents = clusters2process)

# Run the standard workflow for visualization and clustering ------------
new.object <- FindVariableFeatures(object = new.object, selection.method = "vst", nfeatures = 2000)
new.object <- ScaleData(new.object, verbose = F) 
new.object <- RunPCA(new.object, npcs = 30, verbose = FALSE)
new.object <- RunUMAP(new.object, reduction = "pca", dims = 1:20)
new.object <- FindNeighbors(new.object, reduction = "pca", dims = 1:20)
new.object <- FindClusters(new.object, resolution = 0.5)

# get meta data from the integrated object --------------------------------
meta_tab <- new.object@meta.data
meta_tab$barcode <- rownames(meta_tab)

# get DEGs and markers ----------------------------------------------------------------
new.object@misc[["degs"]] <- FindAllMarkers(object = new.object, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
new.object@misc[["markers"]] <- FindAllMarkers(object = new.object, test.use = "roc", only.pos = T, return.thresh = 0.5)
saveRDS(object = new.object, file = paste0(dir_out, case_id_tmp, ".Tummor_Normal.Integrated.20200117.v1.", cell_type_group_tmp, ".StandardReclustered.", run_id, ".RDS"))

# plot dimention reduction ---------------------------------------------
file2write <- paste(dir_out, case_id_tmp, ".Tummor_Normal.Integrated.20200117.v1.", cell_type_group_tmp, ".StandardReclustered.", run_id, ".DimPlot.ByCluster.png", sep="")
png(file = file2write, width = 850, height = 800, res = 150)
p <- DimPlot(new.object, reduction = "umap", group.by = "seurat_clusters", label = T)
p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank())

print(p)
dev.off()

# plot dimention reduction by sample---------------------------------------------
file2write <- paste(dir_out, case_id_tmp, ".Tummor_Normal.Integrated.20200117.v1.", cell_type_group_tmp, ".StandardReclustered.", run_id, ".DimPlot.BySample.png", sep="")
png(file = file2write, width = 850, height = 800, res = 150)
p <- DimPlot(new.object, reduction = "umap", group.by = "orig.ident", label = T)
p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank())

print(p)
dev.off()


# plot dimention reduction by sample---------------------------------------------
file2write <- paste(dir_out, case_id_tmp, ".Tummor_Normal.Integrated.20200117.v1.", cell_type_group_tmp, ".StandardReclustered.", run_id, ".DimPlot.ByClusterSample.png", sep="")
png(file = file2write, width = 2300, height = 800, res = 150)
p <- DimPlot(new.object, reduction = "umap", group.by = "seurat_clusters", label = T, split.by = "orig.ident")
p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank())

print(p)
dev.off()

