# Yige Wu @WashU Oct 2019
## for conducting pathway analysis between tumor subclusters


# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
packages = c(
  "clusterProfiler",
  "ReactomePA",
  "org.Hs.eg.db"
)

for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    BiocManager::install(pkgs = pkg_name_tmp, update = F)
  }
  library(package = pkg_name_tmp, character.only = T)
}

# set run id --------------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)


# set output directory ----------------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# input DEG table ---------------------------------------------------------
deg_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/process_nonimmune_cells_on_cluster/20191022.v1/NonImmune.DEGs.Pos.txt", data.table = F)


# input cluster 2 cell type table -----------------------------------------
cluster2celltype_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Type_Assignment/ccRCC_snRNA_Downstream_Processing - NonImmuneCluster2CellType.20191022.v1.tsv", data.table = F)

# run REACTOME PA ---------------------------------------------------------
for (cluster_tmp in cluster2celltype_tab$Cluster) {
  degs_cluster_tmp <- deg_tab$gene[deg_tab$cluster == cluster_tmp]
  degs_cluster_tmp
  eg = bitr(degs_cluster_tmp, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
  
  x <- enrichPathway(gene=eg$ENTREZID, pvalueCutoff=0.05, readable=T)
  write.table(x = as.data.frame(x), file = paste0(dir_out, "Cluster", cluster_tmp, "_", cluster2celltype_tab$Enriched_Cell_Type_Abbr[cluster2celltype_tab$Cluster == cluster_tmp], "_DEG_Enriched_Pathways.tsv"), quote = F, row.names = F, sep = "\t")
  
  file2write <- paste0(dir_out, "Cluster", cluster_tmp, "_", cluster2celltype_tab$Enriched_Cell_Type_Abbr[cluster2celltype_tab$Cluster == cluster_tmp], "_DEG_Enriched_Pathways.", run_id, ".png")
  png(file = file2write, width = 1500, height = 600, res = 150)
  p <- dotplot(x, showCategory=15)
  print(p)
  dev.off()
}




