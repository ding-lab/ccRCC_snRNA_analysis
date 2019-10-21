# Yige Wu @WashU Sep 2019
## make scatterplot to compare xCell cell type quantification with snRNA-based clutering

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/ccRCC_single_cell_shared.R")


# set parameters ----------------------------------------------------------
version_tmp <- 1
seurat_clusters2sum <- list()
seurat_clusters2sum[["MacroMono"]] <- c(0, 2, 6, 10, 13, 14)
seurat_clusters2sum[["Tcell"]] <- c(8)
# seurat_clusters2sum[["pDC"]] <- c(15)
seurat_clusters2sum[["Endothelial"]] <- c(11, 12)
seurat_clusters2sum[["Epithelial"]] <- c(1, 3, 4, 5, 7, 9)

# input meta data table ---------------------------------------------------
meta_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/sample_info/make_meta_data/meta_data.20190924.v1.tsv", data.table = F)

# input xCell summarized result by major cell clusters in snRNA data --------------------
xcell_value_sum_tab <- fread("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/sample_info/summarize_xCell/xCell_value_sum.20190924.v1.tsv")
xcell_value_sum_tab2merge <- melt(xcell_value_sum_tab, id.vars = c("SampID"))
colnames(xcell_value_sum_tab2merge) <- c("Specimen.ID.bulk", "Cell_Group", "xCell_value")

# Input snRNA barcodes per cluster ----------------------------------------
sn_cell_num_tab <-  fread("./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/analysis_results/scRNA/intergration/integrate_seurat_objects_20190914_v1/snRNA_bc_num_per_cluster.20190924.v1.tsv")
sn_cell_num_tab$Cell_Group <- ""

# get barcode number sum by sample ----------------------------------------
sn_cell_num_per_specimen <- sn_cell_num_tab %>%
  group_by(SampID) %>%
  summarize(Barcode_Sum_Specimen = sum(Barcode_Num))
sn_cell_num_per_specimen
  
# summarize snRNA barcodes by major cell clusters---------------------------------------------------------------
sn_cell_sum_tab <- NULL
for (cell_cluster_tmp in names(seurat_clusters2sum)) {
  sn_cell_tmp_df <- sn_cell_num_tab %>%
    filter(seurat_clusters %in% seurat_clusters2sum[[cell_cluster_tmp]]) %>%
    group_by(SampID) %>%
    summarize(Barcode_Num = sum(Barcode_Num)) %>%
    mutate(seurat_clusters = cell_cluster_tmp)
  sn_cell_tmp_df <- sn_cell_tmp_df %>%
    select(SampID, seurat_clusters, Barcode_Num) %>%
    mutate(Cell_Group = cell_cluster_tmp)
  sn_cell_sum_tab <- rbind(sn_cell_sum_tab, sn_cell_tmp_df)
  
  sn_cell_num_tab$Cell_Group[sn_cell_num_tab$seurat_clusters %in% seurat_clusters2sum[[cell_cluster_tmp]]] <- cell_cluster_tmp
}

sn_perc_tab <- rbind(sn_cell_sum_tab, sn_cell_num_tab)
sn_perc_tab <- merge(sn_perc_tab, sn_cell_num_per_specimen, by = c("SampID"))
sn_perc_tab <- sn_perc_tab %>%
  mutate(Barcode_Perc = Barcode_Num/Barcode_Sum_Specimen) %>%
  select(SampID, Barcode_Perc, seurat_clusters, Cell_Group)

# merge xCell with snRNA barcode number -----------------------------------
sup_tab <- meta_tab
sup_tab <- merge(sup_tab, sn_perc_tab, by.x = c("Specimen.ID.snRNA"), by.y = c("SampID"), all.y = T)
sup_tab <- merge(sup_tab, xcell_value_sum_tab2merge, by = c("Specimen.ID.bulk", "Cell_Group"), all.x = T)
sup_tab <- sup_tab %>%
  mutate(cluter_within_cell_group = !(seurat_clusters == Cell_Group))

current.cluster.ids <- c('0', "2", '6',
                         '8', '9', '10',
                         '13',"14", "15",
                         "1",
                         "3",
                         "4",
                         "5",
                         "7",
                         "11",
                         "12",
                         "Endothelial",
                         "Epithelial",
                         "MacroMono",
                         "Tcell")
new.cluster.ids<- c('MRC1_MacroMono_C0', 
                    "ITGAX_MacroMono_C2", 
                    'MacroMono_C6',
                    'Tcell_C8', 
                    'CDH2_ccRCC_C9', 
                    'PDL1_PDL2_Macrophage_C10',
                    'Proliferating_MacroMono_C13', ## PTPRC/CD45, MRC1/CD206+, CSF1R/CD115+
                    "Macrophage_C14", 
                    "PDL1_PDL2_pDC_C15",
                    "MET_ccRCC_C1",
                    "CA9_ccRCC_C3",
                    "CA9_MET_ccRCC_C4",
                    "HNF1B_ENPP3_ccRCC_C5",
                    "Mesenchymal_C7",
                    "Endothelial_Myofibroblast_C11",
                    "Endothelial_Podocyte_C12",
                    "C11+C12",
                    "C1+C3+C4+C5+C7+C9",
                    "C0+C2+C6+C10+C13+C14",
                    "C8")
sup_tab$cluster_text <- plyr::mapvalues(sup_tab$seurat_clusters, from= current.cluster.ids,to = new.cluster.ids)
write.table(x = sup_tab, file = paste0(makeOutDir(), "snRNA_cell_type_perc_w_xCell", ".", format(Sys.Date(), "%Y%m%d") , ".v", version_tmp , ".tsv"), quote = F, row.names = F, sep = "\t")

# make scatterplot --------------------------------------------------------
## x axis: seurat clusters + sum
## y axis: xCell summarized result by major cell clusters
library(ggrepel)
version_tmp <- 1

tab2p <- sup_tab %>%
  filter(!is.na(xCell_value))

tab2p_text <- tab2p %>%
  filter(Specimen.ID.snRNA == unique(tab2p$Specimen.ID.snRNA)[1])

p <- ggplot()
p <- p + geom_point(data = tab2p, mapping = aes(x = Barcode_Perc, y = xCell_value, color = Case.ID))
p <- p + geom_line(data = tab2p, mapping = aes(x = Barcode_Perc, y = xCell_value, group = seurat_clusters))
p <- p + geom_text_repel(data = tab2p_text, 
                                  mapping = aes(x = Barcode_Perc, 
                                                y = xCell_value, 
                                                label = cluster_text), 
                         force = 4, segment.size = 0.3, segment.alpha = 0.8, angle = 90, size = 3, direction = "x", nudge_y = 0.5, segment.color = "green")
p <- p + facet_grid(cluter_within_cell_group~Cell_Group, space = "fixed", scales = "free")
pdf(file = paste0(makeOutDir(), "scatterplot_Barcode_Num_vs_xCell", ".", format(Sys.Date(), "%Y%m%d") , ".v", version_tmp , ".pdf"), width = 12, height = 6, useDingbats = F)
print(p)
dev.off()


