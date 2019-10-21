# Yige Wu @WashU Sep 2019
## for isolating the non-immune cell clusters and re-do clustering
## plotting workflow: 1) use oberservatoins.txt from infercnv and take the average cnv for each chromosome 2) take barcodes information from oberservatoins.txt 3) umap coordinates can be extracted from Seurat object 4) ggplot do mapping
## ref: https://hbctraining.github.io/scRNA-seq/lessons/06_SC_clustering_quality_control.html

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")


# functions ---------------------------------------------------------------
map_infercnv_category <- function(copy_ratio) {
  cnv_cat <- vector(mode = "character", length = length(copy_ratio))
  cnv_cat[is.na(copy_ratio)] <- "Not Available"
  cnv_cat[!is.na(copy_ratio)] <- "Neutral CN"
  cnv_cat[log2(copy_ratio) < -1 ] <- "Deletion"
  cnv_cat[log2(copy_ratio) >= -1 & log2(copy_ratio) < -0.4] <- "Loss"
  cnv_cat[log2(copy_ratio) >= 0.3 & log2(copy_ratio) < 0.7] <- "Gain"
  cnv_cat[log2(copy_ratio) > 0.7] <- "Amplification"
  return(cnv_cat)
}


# set parameters ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)
sample_id <- c("CPT0086350004")
genes2plot <- c("PBRM1", "BAP1", "SETD2", "HIF1A", "CA9", "JAK2")

# input seurat object -----------------------------------------------------
# seurat_object <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/integration/integrate_seurat_objects/20191002.v1/renal_integrated20191002.v1.RDS")
seurat_object <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/integration/integrate_seurat_objects/20190927.v1/renal_integrated.20190927.v1.RDS")

# extract the UMAP coordinates for each cell ------------------------------
umap_tab <- FetchData(seurat_object, vars = c("orig.ident", "ident", "UMAP_1", "UMAP_2"))
umap_tab$barcode_int <- rownames(umap_tab)
umap_tab <- umap_tab %>%
  filter(orig.ident %in% sample_id) %>%
  mutate(barcode = str_split_fixed(string = barcode_int, pattern = "_", n = 2)[,1])

# input infercnv observations ---------------------------------------------
infercnv_observe_mat <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/Resources/snRNA_Processed_Data/inferCNV/outputs/integration5.20191003.v1/CPT0086350004/infercnv.observations.txt", data.table = F)
dim(infercnv_observe_mat)
infercnv_observe_tab <- melt(infercnv_observe_mat, id.vars = c("V1"))
summary(log2(infercnv_observe_tab$value))

infercnv_observe_tab <- infercnv_observe_tab %>%
  rename(gene_symbol = V1) %>%
  filter(gene_symbol %in% genes2plot) %>%
  mutate(barcode = str_split_fixed(string = variable, pattern = "_", n = 2)[,1]) %>%
  rename(copy_ratio = value) %>%
  select(gene_symbol, barcode, copy_ratio)

for (gene_tmp in genes2plot) {
  
}
# add copy number info to each cell ---------------------------------------
gene_tmp <- "HIF1A"
infercnv_observe_gene_tab <- infercnv_observe_tab %>%
  filter(gene_symbol == gene_tmp) 

tab2p <- umap_tab
tab2p <- merge(tab2p, infercnv_observe_gene_tab, by = c("barcode"), all.x = T)
tab2p$cnv_cat <- map_infercnv_category(copy_ratio = tab2p$copy_ratio)
tab2p$cnv_cat %>% table()
summary(log2(tab2p$copy_ratio))

# plot a UMAP plot for copy number metric ---------------------------------
ggplot(tab2p,
       aes(UMAP_1, UMAP_2)) +
  geom_point(aes(color=cnv_cat), 
             alpha = 0.7) +
  scale_fill_hue()

+
  geom_text(data=umap_label, 
            aes(label=ident, x, y)) +
  ggtitle(qc)

map(metrics, function(qc){

}) %>%
  plot_grid(plotlist = .)
