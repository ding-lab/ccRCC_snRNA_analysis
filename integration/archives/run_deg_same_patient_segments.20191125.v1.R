# Yige Wu @WashU Nov 2019
## for running DEG analysis on the integrated object for multiple pieces per case

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
case_ids <- c("C3N-00733", "C3N-01200")
case_ids <- c("C3N-00733")

case_id_tmp <- "C3N-00733"
for (case_id_tmp in case_ids) {
  ## input seurat object
  if (case_id_tmp == "C3N-00733") {
    seurat_obj_path <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/integrate_same_patient_segments/20191125.v1/", case_id_tmp, ".Tummor_Segments.Integrated.20191125.v1.RDS")
  }
  if (case_id_tmp == "C3N-01200") {
    seurat_obj_path <- paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/integrate_same_patient_segments/20191122.v1/", case_id_tmp, ".Tummor_Segments.Integrated.20191122.v1.RDS")
  }
  seurat_obj <- readRDS(file = seurat_obj_path)
  DefaultAssay(seurat_obj) <- "RNA"
  
  renal.markers <- FindAllMarkers(object = seurat_obj, only.pos = TRUE, min.pct = 0.25, thresh.use = 0.25)
  renal.markers %>%
    colnames()
  renal.markers <- renal.markers[, c("gene", "cluster", "p_val_adj", "p_val", "avg_logFC", "pct.1", "pct.2")]
  write.table(renal.markers, file = paste0(dir_out, case_id_tmp, ".Tumor_Segments.DEGs.Pos.txt"), quote = F, sep = "\t", row.names = F)
  
  
  list_DEGs_by_cluster <- list()
  for (i in unique(renal.markers$cluster)) {
    df2write <- renal.markers %>%
      filter(cluster == i) %>%
      arrange(desc(avg_logFC))
    list_DEGs_by_cluster[[i]] <- df2write
  }
  file2write <- paste0(dir_out, case_id_tmp, ".Tumor_Segments.DEGs.Pos.", run_id, ".xlsx")
  write.xlsx(list_DEGs_by_cluster, file = file2write)
}