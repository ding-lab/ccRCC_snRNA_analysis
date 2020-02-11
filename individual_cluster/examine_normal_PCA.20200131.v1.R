# Yige Wu @WashU Jan. 2020
## For the cell type assignment of the normal sample


# source ------------------------------------------------------------------
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")

# set run id  ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)

# set aliquot id ----------------------------------------------------------
snRNA_aliquot_id_tmp <- "CPT0075170013"

# input seurat processing summary ------------------------------------------------
seurat_summary <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/summary/ccRCC_snRNA_Downstream_Processing - Seurat_Preprocessing.20200116.v1.tsv", data.table = F)
seurat_summary2process <- seurat_summary %>%
  filter(Cellranger_reference_type == "pre-mRNA") %>%
  filter(Proceed_for_downstream == "Yes") %>%
  filter(FACS == "") %>%
  filter(Aliquot %in% snRNA_aliquot_id_tmp) %>%
  mutate(Path_seurat_object = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                     "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                     "/", Aliquot, FACS, "_processed.rds")) %>%
  mutate(Path_deg_table = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                 "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                 "/", Aliquot, FACS, ".DEGs.Pos.txt"))

# input marker table ------------------------------------------------------
gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Kidney_Markers/RCC_marker_gene_table_and_literature_review - Gene2CellType_Tab_w.HumanProteinAtlas.20200131.v1.tsv")

# set output directory ----------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input seurat object -----------------------------------------------------
seurat_obj_path <- seurat_summary2process$Path_seurat_object[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
seurat_obj <- readRDS(file = seurat_obj_path)

# Re-do the scaling to include the cell type marker genes -----------------
all.genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, verbose = F, features = all.genes)

# run the PCA
genes_for_pca <- unique(c(seurat_obj@assays$RNA@var.features, gene2cellType_tab$Gene))
length(genes_for_pca) # 2140
seurat_obj <- RunPCA(seurat_obj, do.print = F, features = genes_for_pca)
# Warning messages:
#   1: In PrepDR(object = object, features = features, verbose = verbose) :
#   The following 3 features requested have not been scaled (running reduction without them): LINC01272, GPX1, RP11-1143G9.4
# 2: In PrepDR(object = object, features = features, verbose = verbose) :
#   The following 28 features requested have zero variance (running reduction without them): ATP6V1G3, CLDN8, AVPR2, CD79B, TCL1A, VPREB3, AHSP, ALAS2, HBA1, HBA2, HBB, HBD, HBE1, PRDX2, TPSAB1, TPSB2, CXCR3, ENHO, MS4A7, CTSD, CTSG, KLRB1, XCL1, XCL2, GZMK, CD3D, SEMA3G, SFRP2

# cluster
seurat_obj <- FindNeighbors(object = seurat_obj, dims = 30)
seurat_obj <- FindClusters(object = seurat_obj, resolution = 0.5)

# run UMAP
seurat_obj <- RunUMAP(seurat_obj, dims = 30)

# save object so far
saveRDS(seurat_obj,file = paste0(dir_out, snRNA_aliquot_id_tmp, "_processed.rds", sep=""))

p <- DimPlot(object = seurat_obj, label = TRUE) + NoLegend()
p

dim(seurat_obj@assays$RNA@scale.data)
dim(seurat_obj@assays$RNA@data)
head(seurat_obj@assays$RNA@data)

# identify differentially expressed cell type marker genes ----------------
## input DEG
deg_tab_path <- seurat_summary2process$Path_deg_table[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
deg_tab_path
deg_tab <- fread(input = deg_tab_path, data.table = F)

cluster_tmp <- 3
cluster_degs <- deg_tab$gene[deg_tab$cluster == cluster_tmp & deg_tab$p_val_adj < 0.1]
## cannot find much cell type marker genes within the DEGs

gene2cellType_cluster3 <- gene2cellType_tab %>%
  filter(Gene %in% cluster_degs)

# get cluster-defining genes---------------------------------------


marker_genes <- FindAllMarkers(object = seurat_obj, test.use = "roc", only.pos = T, return.thresh = 0.5)

cluster_tmp <- 11
cluster_markers <- marker_genes %>%
  filter(cluster == cluster_tmp) %>%
  filter(gene %in% gene2cellType_tab$Gene)

## write as excel table
list_DEGs_by_cluster <- list()
# list_DEGs_by_cluster[["README"]] <- cluster2celltype_tab_tmp
for (i in unique(marker_genes$cluster)) {
  df2write <- marker_genes %>%
    filter(cluster == i) %>%
    filter(power > 0)
  df2write$Cell_Type_Group <- mapvalues(x = ifelse(df2write$gene %in% gene2cellType_tab$Gene, df2write$gene, NA), from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type_Group)
  df2write$Cell_Type1 <- mapvalues(x = ifelse(df2write$gene %in% gene2cellType_tab$Gene, df2write$gene, NA), from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type1)
  df2write$Cell_Type2 <- mapvalues(x = ifelse(df2write$gene %in% gene2cellType_tab$Gene, df2write$gene, NA), from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type2)
  df2write$Cell_Type3 <- mapvalues(x = ifelse(df2write$gene %in% gene2cellType_tab$Gene, df2write$gene, NA), from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type3)
  
  list_DEGs_by_cluster[[as.character(i)]] <- df2write
}
file2write <- paste0(dir_out, snRNA_aliquot_id_tmp, ".AllCluster.ROCTestMarkers.Pos.", run_id, ".xlsx")
write.xlsx(list_DEGs_by_cluster, file = file2write)


# Dotplot -----------------------------------------------------------------
genes2plot <- gene2cellType_tab$Gene
length(genes2plot)
# genes2plot <- intersect(rownames(seurat_obj@assays$RNA@scale.data), genes2plot)
# length(genes2plot)
p <- DotPlot(object = seurat_obj, features = genes2plot, col.min = 0)
p$data$gene_cell_type_group <- plyr::mapvalues(p$data$features.plot, from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type_Group)
p$data$gene_cell_type1 <- plyr::mapvalues(p$data$features.plot, from = gene2cellType_tab$Gene, to = gene2cellType_tab$Cell_Type1)
p <- p  + RotatedAxis()
p <- p + facet_grid(.~gene_cell_type_group + gene_cell_type1, scales = "free", space = "free", drop = T)
p <- p + theme(panel.spacing = unit(0, "lines"),
               strip.background = element_blank(),
               panel.border = element_rect(colour = "black"),
               strip.text.x = element_text(angle = 90, vjust = 0.5, size = 10),
               axis.text.x = element_text(size = 9),
               strip.placement = "outside")

file2write <- paste0(dir_out, snRNA_aliquot_id_tmp,".Individual_Clustered.Dotplot.CellTypeMarkers.", run_id, ".png")
png(file = file2write, width = 4000, height = 2000, res = 150)
print(p)
dev.off()


