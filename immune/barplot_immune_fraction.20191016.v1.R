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

# input integrated data ---------------------------------------------------
object2plot <- readRDS("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/analysis_results/integration/integrate_seurat_objects/20191015.v1/Renal_Integrated.20191015.v1.RDS")
DefaultAssay(object2plot) <- "RNA"
aliquot_ids <- unique(object2plot@meta.data$orig.ident)

# input cluster cell type assignment --------------------------------------
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - AllCluster2Cell_Type.20191015.v2.tsv", data.table = F)
immune_clusters <- cluster2celltype_tab %>%
  filter(Is_Immune == "Yes") %>%
  select(Cluster)
immune_clusters <- immune_clusters$Cluster


# get number of barcodes per cluster --------------------------------------
sn_cell_num_tab <- data.frame(object2plot@meta.data %>%
                                select(orig.ident, seurat_clusters) %>%
                                table())

colnames(sn_cell_num_tab) <- c("snRNA_Aliquot_ID", "Cluster", "Num_Cluster_Barcode")
sn_cell_num_tab <- merge(sn_cell_num_tab, cluster2celltype_tab, by = c("Cluster"), all.x = T)

sn_cell_num_tab <- sn_cell_num_tab %>%
  filter(Is_Immune == "Yes")

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





