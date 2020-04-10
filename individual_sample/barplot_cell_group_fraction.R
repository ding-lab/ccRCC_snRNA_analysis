# Yige Wu @WashU Feb 2020
## for plotting the fraction of immune/stroma/nephron_epithelium

# source ------------------------------------------------------------------
setwd(dir = "~/Box/")
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")

# set run id  ----------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)


# set output directory ----------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input seurat processing summary ------------------------------------------------
## the following path might have issues with the spaces in the file name
seurat_summary <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/summary/ccRCC_snRNA_Downstream_Processing - Seurat_Preprocessing.20200207.v1.tsv", data.table = F)
seurat_summary2process <- seurat_summary %>%
  filter(Cellranger_reference_type == "pre-mRNA") %>%
  filter(Proceed_for_downstream == "Yes") %>%
  filter(FACS == "") %>%
  mutate(Path_seurat_object = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                     "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                     "/", Aliquot, FACS, "_processed.rds"))

# get number of barcodes per cluster --------------------------------------
## get number of barcode per cluster per sample
barcodes_by_cluster_df <- NULL
for (snRNA_aliquot_id_tmp in unique(seurat_summary2process$Aliquot)) {
  ## input seurat object
  seurat_obj_path <- seurat_summary2process$Path_seurat_object[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
  seurat_obj <- readRDS(file = seurat_obj_path)
  
  ## extract meta data and count # barcodes per cluster
  barcodes_by_cluster_aliquot <- data.frame(seurat_obj@meta.data %>%
                                  select(orig.ident, seurat_clusters) %>%
                                  table())
  ## rename column names
  colnames(barcodes_by_cluster_aliquot) <- c("Aliquot", "Cluster", "Num_Cluster_Barcodes")
  
  ## get the sum of barcodes for this sample
  barcodes_aliquot_sum <- sum(barcodes_by_cluster_aliquot$Num_Cluster_Barcodes)
  barcodes_by_cluster_aliquot$Num_Aliquot_Barcodes <- barcodes_aliquot_sum
  
  ## bind with the super table
  barcodes_by_cluster_df <- rbind(barcodes_by_cluster_aliquot, barcodes_by_cluster_df)
}

## get fraction of barcodes per cluster divided by all barcodes in each sample
barcodes_by_cluster_df$Frac_Cluster_Barcodes <- barcodes_by_cluster_df$Num_Cluster_Barcodes/barcodes_by_cluster_df$Num_Aliquot_Barcodes

## write table
file2write <- paste0(dir_out, "Barcodes_By_Cluster.", run_id, ".tsv")
write.table(x = barcodes_by_cluster_df, file = file2write, quote = F, row.names = F, sep = "\t")

# annotate cluster to cell type group -------------------------------------
## input cluster to cell type table (most updated)
cluster2celltype_df <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - Individual.AllCluster2Cell_Type.20200207.v2.tsv", data.table = F)
## merge by aliquot id and cluster name
barcodes_by_cellgroup_df <- merge(barcodes_by_cluster_df, cluster2celltype_df, by = c("Aliquot", "Cluster"), all.x = T)

## sum the barcode fraction by cell group
frac_barcodes_by_cellgroup <- barcodes_by_cellgroup_df %>%
  group_by(Aliquot, Most_Enriched_Cell_Group) %>%
  summarise(Frac_CellGroup_Barcodes = sum(Frac_Cluster_Barcodes))

# make barplot for general 3 cell groups ------------------------------------------------------------
plot_df <- frac_barcodes_by_cellgroup
## label sample to case
plot_df$Case <- mapvalues(x = plot_df$Aliquot, from = seurat_summary2process$Aliquot, to = seurat_summary2process$Case)
## label sample to tumor and normal
plot_df$Sample_Type <- mapvalues(x = plot_df$Aliquot, from = seurat_summary2process$Aliquot, to = seurat_summary2process$Sample_Type)

## ggplot
p <- ggplot()
p <- p + geom_col(data = plot_df, mapping = aes(x = Aliquot, y = Frac_CellGroup_Barcodes, fill = Most_Enriched_Cell_Group), position = "stack")
p <- p + facet_grid(.~Case + Sample_Type, scales = "free", space = "free", drop = T)
p <- p + theme(panel.spacing = unit(0, "lines"),
               # strip.background = element_blank(),
               axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5),
               strip.text.x = element_text(angle = 90, vjust = 0.5, size = 9, face = "bold"),
               panel.border = element_rect(color = "black", fill = NA, size = 1),
               strip.placement = "outside")
p
file2write <- paste0(dir_out, "Cell_Group_Composition.", run_id, ".png")
png(file = file2write, width = 3000, height = 1000, res = 150)
print(p)
dev.off()

# make barplot for immune cell groups ------------------------------------------------------------
plot_df1 <- frac_barcodes_by_cellgroup %>%
  filter(!(Most_Enriched_Cell_Group %in% c("Immune")))
plot_df1 <- as.data.frame(plot_df1)

plot_df2 <- barcodes_by_cellgroup_df %>%
  filter(Most_Enriched_Cell_Group %in% c("Immune")) %>%
  select(Aliquot, Most_Enriched_Cell_Type1, Frac_Cluster_Barcodes) %>%
  rename(Frac_CellGroup_Barcodes = Frac_Cluster_Barcodes) %>%
  rename(Most_Enriched_Cell_Group = Most_Enriched_Cell_Type1)
plot_df2 <- as.data.frame(plot_df2)

plot_df <- base::rbind(plot_df2, plot_df1)

## label sample to case
plot_df$Case <- mapvalues(x = plot_df$Aliquot, from = seurat_summary2process$Aliquot, to = seurat_summary2process$Case)
## label sample to tumor and normal
plot_df$Sample_Type <- mapvalues(x = plot_df$Aliquot, from = seurat_summary2process$Aliquot, to = seurat_summary2process$Sample_Type)


## ggplot
p <- ggplot()
p <- p + geom_col(data = plot_df, mapping = aes(x = Aliquot, y = Frac_CellGroup_Barcodes, fill = Most_Enriched_Cell_Group), position = "stack")
p <- p + facet_grid(.~Case + Sample_Type, scales = "free", space = "free", drop = T)
p <- p + theme(panel.spacing = unit(0, "lines"),
               # strip.background = element_blank(),
               axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5),
               strip.text.x = element_text(angle = 90, vjust = 0.5, size = 9, face = "bold"),
               panel.border = element_rect(color = "black", fill = NA, size = 1),
               strip.placement = "outside")
p
file2write <- paste0(dir_out, "Immune_Cell_Group_Composition.", run_id, ".png")
png(file = file2write, width = 3000, height = 1000, res = 150)
print(p)
dev.off()
