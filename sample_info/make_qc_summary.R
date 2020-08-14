# Yige Wu @ WashU 2019 Nov
## make summary table listing estimated number of cells, mean reads per cell, median genes per cell and nuber of cells passed QC

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input seurat processing summary ------------------------------------------------
summary_df <- fread(data.table = F, input = "./Resources/snRNA_Processed_Data/scRNA_auto/summary/Seurat_Preprocessing.20200701.v1.tsv")
## input case id and aliquot id
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")

# process summary table ---------------------------------------------------
summary_df <- summary_df %>%
  filter(Proceed_for_downstream == "Yes") %>%
  filter(Cellranger_reference_type == "pre-mRNA") %>%
  filter(FACS == "") %>%
  mutate(Path_metrics_summary = paste0("./Resources/snRNA_Processed_Data/Cell_Ranger/outputs/", 
                                       Aliquot, "/outs/metrics_summary.csv"))

# extract info from metrics_summary.csv -----------------------------------
snRNA_aliquot_id_tmp <- "CPT0001220012"
metrics_summary_ext_tab <- NULL
for (snRNA_aliquot_id_tmp in unique(summary_df$Aliquot)) {
  path_metrics_summary_tmp <- summary_df$Path_metrics_summary[summary_df$Aliquot == snRNA_aliquot_id_tmp]
  metrics_summary_tmp <- fread(input = path_metrics_summary_tmp, data.table = F)
  metrics_summary_ext_tmp <- metrics_summary_tmp %>%
    select(`Estimated Number of Cells`, `Mean Reads per Cell`, `Median Genes per Cell`)
  metrics_summary_ext_tmp$Number_of_Cells_Passed_QC <- summary_df$No.valid_barcode[summary_df$Aliquot == snRNA_aliquot_id_tmp]
  metrics_summary_ext_tmp <- cbind(data.frame(Aliquot = snRNA_aliquot_id_tmp), metrics_summary_ext_tmp)
  metrics_summary_ext_tab <- rbind(metrics_summary_ext_tmp, metrics_summary_ext_tab)
}
metrics_summary_ext_tab$Id_Case <- mapvalues(x = metrics_summary_ext_tab$Aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Case))
metrics_summary_ext_tab$Id_Aliquot_WU <- mapvalues(x = metrics_summary_ext_tab$Aliquot, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))

# write output ------------------------------------------------------------
write.table(x = metrics_summary_ext_tab, file = paste0(dir_out, "QC_Metrics_Summary.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)
