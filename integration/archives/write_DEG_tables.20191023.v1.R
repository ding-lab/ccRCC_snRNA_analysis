# Yige Wu @WashU Oct 2019
## for outputing DEG tables

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")


# set run id --------------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)


# set output directory ----------------------------------------------------------
integration_id <- "20191021.v1"
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input DEG ---------------------------------------------------------------
deg_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/integrate_seurat_objects/20191021.v1/Renal.DEGs.Pos.txt", data.table = F)


# input cluster 2 cell type table -----------------------------------------
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - AllCluster2Cell_Type.20191022.v1.tsv", data.table = F)

# write tables ------------------------------------------------------------
list_DEGs_by_cluster <- list()
list_DEGs_by_cluster[["README"]] <- cluster2celltype_tab
for (i in unique(deg_tab$cluster)) {
  df2write <- deg_tab %>%
    filter(cluster == i) %>%
    arrange(desc(avg_logFC))
  list_DEGs_by_cluster[[as.character(i)]] <- df2write
}
file2write <- paste0(dir_out, "Renal.DEGs.Pos.", run_id, ".xlsx")
write.xlsx(list_DEGs_by_cluster, file = file2write)
