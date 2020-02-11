# Yige Wu @WashU Sep 2019
## for integrating  snRNA datasets belong to the same patient

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

# set case ids to be processed --------------------------------------------
case_id <- c("C3N-01200")

# input seurat processing summary ------------------------------------------------
seurat_summary <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/summary/ccRCC_snRNA_Downstream_Processing - Seurat_Preprocessing.20200204.v1.tsv", data.table = F)
seurat_summary2process <- seurat_summary %>%
  filter(Proceed_for_downstream == "Yes") %>%
  filter(Cellranger_reference_type == "pre-mRNA") %>%
  filter(Sample_Type == "Tumor") %>%
  filter(Case %in% case_id) %>%
  filter(FACS == "") %>%
  mutate(Path_seurat_object = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                     "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                     "/", Aliquot, FACS, "_processed.rds"))
seurat_summary2process$Path_seurat_object

# set aliquot ids to be processed -----------------------------------------
snRNA_aliquot_ids <- seurat_summary2process$Aliquot

# Input seurat objects -----------------------------------------------------
renal.list <- list()
for (i in 1:nrow(seurat_summary2process)) {
  sample_id_tmp <- seurat_summary2process$Aliquot[i]
  seurat_obj_path <- seurat_summary2process$Path_seurat_object[i]
  seurat_obj <- readRDS(file = seurat_obj_path)
  
  seurat_obj$orig.ident  <- sample_id_tmp
  DefaultAssay(seurat_obj) <- "RNA" 
  seurat_obj@assays$SCT <- NULL
  seurat_obj@graphs <- list()
  seurat_obj@neighbors <- list()
  seurat_obj@reductions <- list()
  for (col_name_tmp in c("nCount_SCT", "nFeature_SCT", "SCT_snn_res.0.5", "seurat_clusters")) {
    if (col_name_tmp %in% names(seurat_obj@meta.data)) {
      seurat_obj@meta.data[[col_name_tmp]] <- NULL
    }
  }
  renal.list[[sample_id_tmp]] <- seurat_obj
  
}
length(renal.list)
rm(seurat_obj)

#  split the combined object into a list, with each dataset as an element ----------------------------------------
renal.list <- lapply(X = renal.list, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})


# Merge Data -----------------------------------------------------------
renal.merged <- merge(renal.list[[1]], y = renal.list[[2]], project = "Merged")

# Run the standard workflow for visualization and clustering ------------
# The variable features of this assay are automatically set during IntegrateData
renal.merged <- SCTransform(renal.merged, vars.to.regress = c("nCount_RNA","percent.mito"))
renal.merged <- RunPCA(renal.merged, npcs = 30, verbose = FALSE)
renal.merged <- RunUMAP(renal.merged, reduction = "pca", dims = 1:30)
renal.merged <- FindNeighbors(renal.merged, reduction = "pca", dims = 1:30)
renal.merged <- FindClusters(renal.merged, resolution = 0.5)
saveRDS(object = renal.merged, file = paste0(dir_out, case_id, ".Tumor_Segments.Merged.", run_id, ".RDS"))

# Dimplot by sample -------------------------------------------------------
p <- DimPlot(renal.merged, reduction = "umap", group.by = "orig.ident", label = F, label.size	= 5)
p <- p + ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                        axis.text.y=element_blank(),axis.ticks=element_blank(),
                        axis.title.x=element_blank(),
                        axis.title.y=element_blank(),
                        legend.position = "top")
file2write <- paste(dir_out, case_id, ".Tumor_Segments.Merged.Dimplot.BySample.", run_id, ".png", sep="")
png(file = file2write, width = 1200, height = 1300, res = 150)
print(p)
dev.off()

# Identify genes defining each cluster using ROC test------------------------------------
## input cell type marker table
gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/scripts/src/Marker_Gene_Tables/Gene2CellType_Tab.20200203.v1.tsv", data.table = F)
## run FindAllMarkers
markers_roc <- FindAllMarkers(object = renal.merged, test.use = "roc", only.pos = T, return.thresh = 0.5)
## annotate cell type marker genes and write as tsv file
markers_roc$Cell_Type_Group <- mapvalues(x = ifelse(markers_roc$gene %in% gene2cellType_tab$Gene, markers_roc$gene, NA), from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type_Group)
markers_roc$Cell_Type1 <- mapvalues(x = ifelse(markers_roc$gene %in% gene2cellType_tab$Gene, markers_roc$gene, NA), from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type1)
markers_roc$Cell_Type2 <- mapvalues(x = ifelse(markers_roc$gene %in% gene2cellType_tab$Gene, markers_roc$gene, NA), from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type2)
markers_roc$Cell_Type3 <- mapvalues(x = ifelse(markers_roc$gene %in% gene2cellType_tab$Gene, markers_roc$gene, NA), from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type3)
file2write <- paste0(dir_out, case_id, ".Tumor_Segments.Merged.FindAllMarkers.ROC.Pos.", run_id, ".tsv")
write.table(markers_roc, file = file2write, quote = F, sep = "\t", row.names = F)
## write as excel file
list_DEGs_by_cluster <- list()
for (i in unique(markers_roc$cluster)) {
  df2write <- markers_roc %>%
    filter(cluster == i)
  list_DEGs_by_cluster[[as.character(i)]] <- df2write
}
file2write <- paste0(dir_out, case_id, ".Tumor_Segments.Merged.FindAllMarkers.ROC.Pos.", run_id, ".xlsx")
write.xlsx(list_DEGs_by_cluster, file = file2write)

# Identify differentially expressed genes for each cluster using wilcox test------------------------------------
## input cell type marker table
gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/scripts/src/Marker_Gene_Tables/Gene2CellType_Tab.20200203.v1.tsv", data.table = F)
## run FindAllMarkers
markers_wilcox <- FindAllMarkers(object = renal.merged, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
## annotate cell type marker genes and write as tsv file
markers_wilcox$Cell_Type_Group <- mapvalues(x = ifelse(markers_wilcox$gene %in% gene2cellType_tab$Gene, markers_wilcox$gene, NA), from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type_Group)
markers_wilcox$Cell_Type1 <- mapvalues(x = ifelse(markers_wilcox$gene %in% gene2cellType_tab$Gene, markers_wilcox$gene, NA), from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type1)
markers_wilcox$Cell_Type2 <- mapvalues(x = ifelse(markers_wilcox$gene %in% gene2cellType_tab$Gene, markers_wilcox$gene, NA), from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type2)
markers_wilcox$Cell_Type3 <- mapvalues(x = ifelse(markers_wilcox$gene %in% gene2cellType_tab$Gene, markers_wilcox$gene, NA), from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type3)
file2write <- paste0(dir_out, case_id, ".Tumor_Segments.Merged.FindAllMarkers.Wilcox.Pos.", run_id, ".tsv")
write.table(markers_wilcox, file = file2write, quote = F, sep = "\t", row.names = F)
## write as excel file
list_DEGs_by_cluster <- list()
for (i in unique(markers_wilcox$cluster)) {
  df2write <- markers_wilcox %>%
    filter(cluster == i)
  list_DEGs_by_cluster[[as.character(i)]] <- df2write
}
file2write <- paste0(dir_out, case_id, ".Tumor_Segments.Merged.FindAllMarkers.Wilcox.Pos.", run_id, ".xlsx")
write.xlsx(list_DEGs_by_cluster, file = file2write)