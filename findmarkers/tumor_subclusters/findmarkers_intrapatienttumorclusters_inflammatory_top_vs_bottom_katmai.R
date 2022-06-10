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
path_rds <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/Analysis_Results/merging/merge_34_ccRCC_samples/20211005.v1//ccRCC.34samples.Merged.20211005.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")
## input the barcode-manualsubcluster info
barcode2subclusterid_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_tumorsubclusterid/map_barcode_with_manual_tumorsubcluster_id/20210805.v1/Barcode2TumorSubclusterId.20210805.v1.tsv", data.table = F)
enrich_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/assign_tumorcluster_by_msigdb_geneset_scores/20210805.v1/MsigDB_Hallmark.Top15GeneSets.4Module.Enrichment.tsv")

# set parameters ----------------------------------------------------------
## set min.pct
logfc.threshold.run <- 0
min.pct.run <- 0.1
min.diff.pct.run <- 0
## spcify assay
assay_process <- "RNA"
cat(paste0("Assay: ", assay_process, "\n"))
cat("###########################################\n")

# preprocess --------------------------------------------------------------
enrich_df <- enrich_df %>%
  mutate(cluster_name.formatted = gsub(x = cluster_name, pattern = "\\.", replacement = "-"))
cluster_group1_process <- enrich_df$cluster_name.formatted[enrich_df$INFLAMMATORY_RESPONSE_Score >= quantile(x = enrich_df$INFLAMMATORY_RESPONSE_Score, probs = 0.9)]
cluster_group2_process <- enrich_df$cluster_name.formatted[enrich_df$INFLAMMATORY_RESPONSE_Score <= quantile(x = enrich_df$INFLAMMATORY_RESPONSE_Score, probs = 0.1)]

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
srat@meta.data$group_findmarkers[srat@meta.data$group_findmarkers == srat@meta.data$id_aliquot_barcode] <- "others"
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
file2write <- paste0(dir_out, "Inflammatory_score_top_vs_bottom_tumorclusters.", 
                     ".logfc.threshold", logfc.threshold.run, 
                     ".min.pct", min.pct.run,
                     ".min.diff.pct", min.diff.pct.run,
                     ".Assay", assay_process,
                     ".tsv")
write.table(x = deg_df, file = file2write, sep = "\t", quote = F, row.names = F)
cat("finish writing the result!\n\n")


