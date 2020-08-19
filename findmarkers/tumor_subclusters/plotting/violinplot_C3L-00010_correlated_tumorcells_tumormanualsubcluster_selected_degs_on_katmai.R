# Yige Wu @WashU March 2020
## for each individual sample tumor cell reclustered, plot UMAP

# set up libraries and output directory -----------------------------------
## getting the path to the current script
thisFile <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
path_this_script <- thisFile()
## set working directory
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input barcode-to-manual grouping 
barcode2manualsubcluster_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_corrected_by_individual_sample_inspection/20200811.v1/31Aliquot.Barcode2CellType.20200811.v1.tsv", data.table = F)
## load meta data
idmetadata_df <- fread(input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200505.v1/meta_data.20200505.v1.tsv", data.table = F)
## input srat object
srat <- readRDS(file = "./Resources/Analysis_Results/integration/merge_different_cases/merge_tumorcells_C3L-00010_C3L-00416_C3L-00583_C3N-00242_on_katmai/20200817.v1/C3L-00416_C3L-00583_C3L-00010_C3N-00242.Tumor_Segments.Merged.20200817.v1.RDS")
## set genes to plot
genes2plot <- c("CA9", "PPARA", "RORA")

# edit meta data ------------------------------------
## add case id
barcode2manualsubcluster_df$id_case <- mapvalues(x = barcode2manualsubcluster_df$orig.ident, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Case))
barcode2manualsubcluster_df$id_aliquot_wu <- mapvalues(x = barcode2manualsubcluster_df$orig.ident, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
## merge info
metadata_df <- srat@meta.data
metadata_df$barcode_integrated_case <- rownames(srat@meta.data)
cat("Added barcode info from meta data!\n")
metadata_df <- metadata_df %>%
  dplyr::mutate(barcode_individual = str_split_fixed(string = barcode_integrated_case, pattern = "_", n = 3)[,1])
barcode2manualsubcluster_df <- merge(barcode2manualsubcluster_df, metadata_df, 
                                             by.x = c("orig.ident", "individual_barcode"), by.y = c("orig.ident", "barcode_individual"), all.y = T)

## filter down to only tumor cells with manual cluster assigned
barcode2manualsubcluster_df <- barcode2manualsubcluster_df %>%
  dplyr::mutate(Name_Cluster = paste0(id_aliquot_wu, "_C", (Id_TumorManualCluster+1)))
### subset
srat_plot <- subset(srat, cells = barcode2manualsubcluster_df$barcode_integrated_case[!is.na(barcode_integrated_case$Name_Cluster)])
## change meta data
srat_plot@meta.data$Name_Cluster <- mapvalues(x = rownames(metadata_df), from = barcode2manualsubcluster_df$barcode_integrated_case, to = as.vector(barcode2manualsubcluster_df$Name_Cluster))
### set the identities to cluster in the meta data
Idents(object = srat_plot) <- "Name_Cluster"


# plot --------------------------------------------------------------------
for (gene_tmp in genes2plot) {
  p <- VlnPlot(object = srat_plot, features = gene_tmp, group.by = "Name_Cluster", pt.size = 0, ncol = 5)
  p <- p + theme(axis.text.x = element_text(angle = 90))
  p <- p + theme(legend.position = "none")
  png(filename = paste0(dir_out, "C3L-00010_correlated_tumorcells", ".", gene_tmp, ".png"),width = 1000, height = 800, res = 150)
  print(p)
  dev.off()
}
cat("Finished all\n")
