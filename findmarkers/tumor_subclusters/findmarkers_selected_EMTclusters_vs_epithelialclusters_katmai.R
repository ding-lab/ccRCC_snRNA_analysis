# Yige Wu @WashU Apr 2021

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
library(future)
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 7 * 1024^3) # for 7 Gb RAM

# input dependencies ------------------------------------------------------
## input the integrated data
path_rds <- "./Resources/Analysis_Results/merging/merge_35_samples/20210802.v1/RCC.35samples.Merged.20210802.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")
## input the barcode-manualsubcluster info
barcode2subclusterid_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_tumorsubclusterid/map_barcode_with_manual_tumorsubcluster_id/20210805.v1/Barcode2TumorSubclusterId.20210805.v1.tsv", data.table = F)

# set parameters ----------------------------------------------------------
## set min.pct
logfc.threshold.run <- 0
min.pct.run <- 0.1
min.diff.pct.run <- 0.1
## spcify assay
assay_process <- "RNA"
cat(paste0("Assay: ", assay_process, "\n"))
cat("###########################################\n")
## tumor-cell cluster enrichment modules
cluster_group1_process <- c("C3L-00079-T1_C4", "C3L-01302-T1_C1")
cluster_group2_process <- c("C3L-00088-T2_C1", "C3N-00733-T1_C1", "C3L-00416-T2_C1", "C3L-00010-T1_C1", "C3L-00088-T1_C1")

# preprocess --------------------------------------------------------------
## process the barcode info
barcode2subclusterid_df <- barcode2subclusterid_df %>%
  mutate(id_aliquot_barcode = paste0(orig.ident, "_", barcode)) %>%
  mutate(group_findmarkers = ifelse(Cluster_Name %in% cluster_group1_process, "group1",
                                    ifelse(Cluster_Name %in% cluster_group2_process, "group2", "others")))
barcode2subclusterid_df %>%
  select(Cluster_Name, group_findmarkers) %>%
  unique() %>%
  arrange(group_findmarkers) %>%
  select(group_findmarkers) %>%
  table()
table(barcode2subclusterid_df$group_findmarkers)

## process the seurat object
## get original barcode
BC <- srat@meta.data %>% rownames
srat@meta.data$original_barcode <- BC %>% strsplit("_") %>% lapply("[[",1) %>% unlist
head(srat@meta.data$original_barcode)
cat("finish adding the simple barcode!\n")
## make combined id for the seurat meta data
srat@meta.data$id_aliquot_barcode <- paste0(srat@meta.data$orig.ident, "_", srat@meta.data$original_barcode)
head(srat@meta.data$id_aliquot_barcode)
cat("finish adding unique id for each barcode in the seurat object!\n")
srat@meta.data$group_findmarkers <- mapvalues(x = srat@meta.data$id_aliquot_barcode, from = barcode2subclusterid_df$id_aliquot_barcode, to = as.vector(barcode2subclusterid_df$group_findmarkers), warn_missing = F)
table(srat@meta.data$group_findmarkers)

# run findmarkers ---------------------------------------------------------
## run findmarkers
DefaultAssay(srat) <- assay_process
Idents(srat) <- "group_findmarkers" 
deg_df <- FindMarkers(object = srat, test.use = "wilcox", ident.1 = "group1", ident.2 = "group2", only.pos = F,
                      min.pct = min.pct.run, logfc.threshold = logfc.threshold.run, min.diff.pct = min.diff.pct.run, verbose = T)
deg_df$genesymbol_deg <- rownames(deg_df)
deg_df$cellnumber_group1 <- length(which(srat@meta.data$group_findmarkers == "group1"))
deg_df$cellnumber_group2 <- length(which(srat@meta.data$group_findmarkers == "group2"))


# write output ------------------------------------------------------------
## write output
file2write <- paste0(dir_out, "Selected_2EMTclusters_vs_5Epithelialclusters", 
                     ".logfc.threshold", logfc.threshold.run, 
                     ".min.pct", min.pct.run,
                     ".min.diff.pct", min.diff.pct.run,
                     ".Assay", assay_process,
                     ".tsv")
write.table(x = deg_df, file = file2write, sep = "\t", quote = F, row.names = F)
cat("finish writing the result!\n\n")


