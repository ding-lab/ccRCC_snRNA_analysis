# Yige Wu @WashU Oct 2019
## for plotting the fraction of immune cell populations across samples

# source ------------------------------------------------------------------
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")

# set run id  ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)


# set output directory ----------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# set aliquot ids to be processed -----------------------------------------
snRNA_aliquot_ids <- c("CPT0019130004", "CPT0001260013", "CPT0086350004", "CPT0010110013", "CPT0001180011", "CPT0025890002", "CPT0075140002", "CPT0020120013", "CPT0001220012", "CPT0014450005")

# input seurat processing summary ------------------------------------------------
seurat_summary <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/summary/ccRCC_snRNA_Downstream_Processing - Seurat_Preprocessing.20200116.v1.tsv", data.table = F)
seurat_summary2process <- seurat_summary %>%
  filter(Cellranger_reference_type == "pre-mRNA") %>%
  filter(Proceed_for_downstream == "Yes") %>%
  filter(FACS == "") %>%
  filter(Aliquot %in% snRNA_aliquot_ids) %>%
  mutate(Path_seurat_object = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                     "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                     "/", Aliquot, FACS, "_processed.rds"))
  

# input cluster cell type assignment --------------------------------------
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - Individual.AllCluster2Cell_Type.20200116.v1.tsv", data.table = F)

# get number of barcodes per cluster --------------------------------------
sn_cell_sum_tab <- NULL
for (snRNA_aliquot_id_tmp in snRNA_aliquot_ids) {
  ## input seurat object
  seurat_obj_path <- seurat_summary2process$Path_seurat_object[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
  seurat_obj <- readRDS(file = seurat_obj_path)
  
  sn_cell_num_tmp <- data.frame(seurat_obj@meta.data %>%
                                  select(orig.ident, seurat_clusters) %>%
                                  table())
  
  colnames(sn_cell_num_tmp) <- c("snRNA_Aliquot_ID", "Cluster", "Num_Cluster_Barcode")
  
  ## get the cluster2celltype info for current aliquot
  cluster2celltype_tab_tmp <- cluster2celltype_tab %>%
    filter(Aliquot == snRNA_aliquot_id_tmp)
  
  sn_cell_num_tmp <- merge(sn_cell_num_tmp, cluster2celltype_tab_tmp, by = c("Cluster"), all.x = T)
  
  sn_cell_sum_tab <- rbind(sn_cell_num_tmp, sn_cell_sum_tab)
}






sn_cell_sum_tab <- sn_cell_num_tab %>%
  group_by(snRNA_Aliquot_ID) %>%
  summarise(Num_Aliquot_Immune_Cells = sum(Num_Cluster_Barcode))

sn_cell_num_tab <- merge(sn_cell_num_tab, sn_cell_sum_tab, by = c("snRNA_Aliquot_ID"), all.x = T)
sn_cell_num_tab <- sn_cell_num_tab %>%
  mutate(Perc_Cluster_in_Immune = Num_Cluster_Barcode/Num_Aliquot_Immune_Cells)

sn_cell_num_tab %>%
  head()


# plot Perc_Cluster_in_Immune barplot ------------------------------------------------------------
tab2p <- sn_cell_num_tab

p <- ggplot()
p <- p + geom_col(data = tab2p, mapping = aes(x = snRNA_Aliquot_ID, y = Perc_Cluster_in_Immune, fill = Enriched_Cell_Type_Abbr))
p <- p + scale_fill_brewer(palette = "Set1")
p <- p + coord_flip()
p <- p + theme_cowplot()
p1 <- p
file2write <- paste0(dir_out, "Perc_Cluster_in_Immune.", run_id, ".png")
png(file2write, width = 1500, height = 700, res = 150)
print(p)
dev.off()

file2write <- paste0(dir_out, "Perc_Cluster_in_Immune_thin.", run_id, ".png")
png(file2write, width = 1000, height = 500, res = 150)
print(p)
dev.off()

# plot Num_Cluster_Barcode barplot ------------------------------------------------------------
tab2p <- sn_cell_num_tab

p <- ggplot()
p <- p + geom_col(data = tab2p, mapping = aes(x = snRNA_Aliquot_ID, y = Num_Cluster_Barcode, fill = Enriched_Cell_Type_Abbr))
p <- p + scale_fill_brewer(palette = "Set1")
p <- p + coord_flip()
p <- p + scale_y_reverse()
p <- p + theme_cowplot()
p2 <- p
file2write <- paste0(dir_out, "Num_Cluster_Barcode.", run_id, ".png")
png(file2write, width = 1000, height = 500, res = 150)
print(p)
dev.off()

plot_grid(p2, p1)





