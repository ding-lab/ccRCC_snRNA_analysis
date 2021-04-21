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
dir_out1 <- "./Resources/snRNA_Processed_Data/Differentially_Expressed_Genes/Tumor_Subcluster_Group/"
dir.create(dir_out1)
dir_out <- paste0(dir_out1, run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the integrated data
path_rds <- "./Resources/Analysis_Results/merging/RunPCA_UMAP_clustering_32_aliquot/20210318.v1/32_aliquot.Merged.20210318.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")
## input the barcode-manualsubcluster info
barcode2subclusterid_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_tumorsubclusterid/map_barcode_with_manual_tumorsubcluster_id/20201130.v1/Barcode2TumorSubclusterId.20201130.v1.tsv", data.table = F)
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20200902.v1/scrublet.run20200902_adj_cutoff.united_outputs.tsv", data.table = F)
## input enriched gene set module assignment per tumor cluster
enriched_df <- fread(data.table = F, input = "./Resources/Analysis_Results/tumor_subcluster/calculate_scores/assign_tumorcluster_by_msigdb_geneset_scores/20210421.v3/MsigDB_Hallmark.Top15GeneSets.4Module.Enrichment.tsv")

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
modules_process <- c("Cell_cycle", "Immune", "EMT", "mTOR")

# preprocess --------------------------------------------------------------
## revise cluster names
enriched_df <- enriched_df %>%
  mutate(id_manual_cluster = gsub(x = cluster_name, pattern = "\\.", replacement = "-"))
## merge cluster/scrublet/cluster group info
barcode2info_df <- merge(x = barcode2subclusterid_df, y = barcode2scrublet_df, 
                         by.x = c("orig.ident", "barcode", "easy_id"), 
                         by.y = c("Aliquot", "Barcodes", "Aliquot_WU"), all.x = T)
barcode2info_df <- merge(x = barcode2info_df %>%
                           rename(id_manual_cluster = Cluster_Name), 
                         y = enriched_df[, c("id_manual_cluster", modules_process)],
                         by = c("id_manual_cluster"), all.x = T)
## get original barcode
BC <- srat@meta.data %>% rownames
srat@meta.data$original_barcode <- BC %>% strsplit("_") %>% lapply("[[",1) %>% unlist
head(srat@meta.data$original_barcode)
cat("finish adding the simple barcode!\n")
## make combined id for the seurat meta data
srat@meta.data$id_aliquot_barcode <- paste0(srat@meta.data$orig.ident, "_", srat@meta.data$original_barcode)
head(srat@meta.data$id_aliquot_barcode)
cat("finish adding unique id for each barcode in the seurat object!\n")

# for each aliquot, input seurat object and fetch data and write data --------------------
markers_wilcox_df <- NULL
for (module_tmp in modules_process) {
  file2write <- paste0(dir_out,
                       module_tmp, ".TumorManualCluster.DEGs.Wilcox.Minpct", min.pct.run, 
                       ".Logfc", logfc.threshold.run,
                       ".min.diff.pct", min.diff.pct.run,
                       ".tsv")
  if (!file.exists(file2write)) {
    ## make combined id for the barcode2celltype table
    barcode2info_df <- barcode2info_df %>%
      mutate(id_aliquot_barcode = paste0(orig.ident, "_", barcode)) %>%
      mutate(is_module = module_tmp) %>%
      mutate(group_findmarkers = ifelse(!predicted_doublet, ifelse(is_module, "group1", "group2"), "others"))
    
    ## map group label
    srat@meta.data$group_findmarkers <- mapvalues(x = srat@meta.data$id_aliquot_barcode, from = barcode2info_df$id_aliquot_barcode, to = as.vector(barcode2info_df$group_findmarkers), warn_missing = F)
    srat@meta.data$group_findmarkers[srat@meta.data$group_findmarkers == srat@meta.data$id_aliquot_barcode] <- "others"
    table(srat@meta.data$group_findmarkers)
    cat("finish adding group labels\n")
    
    ## run findmarkers
    Idents(srat) <- "group_findmarkers" 
    deg_df <- FindMarkers(object = srat, test.use = "wilcox", ident.1 = "group1", ident.2 = "group2", only.pos = F,
                          min.pct = min.pct.run, logfc.threshold = logfc.threshold.run, min.diff.pct = min.diff.pct.run, verbose = T)
    deg_df$genesymbol_deg <- rownames(deg_df)
    deg_df$module <- module_tmp
    deg_df$cellnumber_moduleenriched <- length(which(srat@meta.data$group_findmarkers == "group1"))
    deg_df$cellnumber_othertumorcells <- length(which(srat@meta.data$group_findmarkers == "group2"))
    
    ## write output
    write.table(x = deg_df, file = file2write, sep = "\t", quote = F, row.names = F)
    cat("finish writing the result!\n\n")
  }
}


