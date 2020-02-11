# Yige Wu @WashU Sep 2019
## for isolating the non-immune cell clusters and re-do clustering
## plotting workflow: 1) use oberservatoins.txt from infercnv and take the average cnv for each chromosome 2) take barcodes information from oberservatoins.txt 3) umap coordinates can be extracted from Seurat object 4) ggplot do mapping
## ref: https://hbctraining.github.io/scRNA-seq/lessons/06_SC_clustering_quality_control.html

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")

# set parameters ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
sample_id <- c("CPT0086350004")
genes2plot <- c("VHL")

# input seurat object -----------------------------------------------------
seurat_object <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/integration/integrate_seurat_objects/20191002.v1/renal_integrated20191002.v1.RDS")

# extract the UMAP coordinates for each cell ------------------------------
umap_tab <- FetchData(seurat_object, vars = c("orig.ident", "ident", "UMAP_1", "UMAP_2"))
umap_tab$barcode <- rownames(umap_tab)
umap_tab <- umap_tab %>%
  filter(orig.ident %in% sample_id)

# input infercnv observations ---------------------------------------------
infercnv_observe_mat <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/snRNA_Processed_Data/inferCNV/outputs/integration5.20191003.v1/CPT0086350004/infercnv.observations.txt", data.table = F)
infercnv_observe_tab <- melt(infercnv_observe_mat, id.vars = c("V1"))
infercnv_observe_tab <- infercnv_observe_tab %>%
  rename(gene_symbol = V1) %>%
  rename(barcode = variable) %>%
  rename(copy_ratio = value) %>%
  filter(gene_symbol %in% genes2plot)

# add copy number info to each cell ---------------------------------------


# plot a UMAP plot for copy number metric ---------------------------------

map(metrics, function(qc){
  ggplot(qc_data,
         aes(UMAP_1, UMAP_2)) +
    geom_point(aes_string(color=qc), 
               alpha = 0.7) +
    scale_color_gradient(guide = FALSE, 
                         low = "grey90", 
                         high = "blue")  +
    geom_text(data=umap_label, 
              aes(label=ident, x, y)) +
    ggtitle(qc)
}) %>%
  plot_grid(plotlist = .)
