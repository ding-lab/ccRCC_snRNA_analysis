# Yige Wu @ WashU 2019 Nov
## make summary table listing estimated number of cells, mean reads per cell, median genes per cell and nuber of cells passed QC

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

# set aliquot ids to be processed -----------------------------------------
snRNA_aliquot_ids <- c("CPT0019130004", "CPT0001260013", "CPT0086350004", "CPT0010110013", "CPT0001180011", "CPT0025890002", "CPT0075140002", "CPT0020120013", "CPT0001220012", "CPT0014450005")

# input seurat processing summary ------------------------------------------------
seurat_summary <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/summary/ccRCC_snRNA_Downstream_Processing - Seurat_Preprocessing.20191021.v1.tsv", data.table = F)
seurat_summary2process <- seurat_summary %>%
  filter(Cellranger_reference_type == "pre-mRNA") %>%
  filter(Proceed_for_downstream == "Yes") %>%
  filter(Aliquot %in% snRNA_aliquot_ids) %>%
  mutate(Path_metrics_summary = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/Cell_Ranger/outputs/", Aliquot, "/outs/metrics_summary.csv"))


# extract info from metrics_summary.csv -----------------------------------
snRNA_aliquot_id_tmp <- "CPT0001220012"
metrics_summary_ext_tab <- NULL
for (snRNA_aliquot_id_tmp in snRNA_aliquot_ids) {
  path_metrics_summary_tmp <- seurat_summary2process$Path_metrics_summary[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
  metrics_summary_tmp <- fread(input = path_metrics_summary_tmp, data.table = F)
  metrics_summary_ext_tmp <- metrics_summary_tmp %>%
    select(`Estimated Number of Cells`, `Mean Reads per Cell`, `Median Genes per Cell`)
  metrics_summary_ext_tmp$Number_of_Cells_Passed_QC <- seurat_summary2process$No.valid_barcode[seurat_summary2process$Aliquot == snRNA_aliquot_id_tmp]
  metrics_summary_ext_tmp <- cbind(data.frame(Aliquot = snRNA_aliquot_id_tmp), metrics_summary_ext_tmp)
  metrics_summary_ext_tab <- rbind(metrics_summary_ext_tmp, metrics_summary_ext_tab)
}
write.table(x = metrics_summary_ext_tab, file = paste0(dir_out, "QC_Metrics_Summary.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)
