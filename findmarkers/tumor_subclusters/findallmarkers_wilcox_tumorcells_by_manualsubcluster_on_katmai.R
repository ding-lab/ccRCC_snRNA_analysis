# Yige Wu @WashU Feb 2020

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
dir_out1 <- "./Resources/snRNA_Processed_Data/Differentially_Expressed_Genes/Tumor_Subclusters/"
dir.create(dir_out1)
dir_out <- paste0(dir_out1, run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the paths for reclustered seurat objects
srat_paths_df <- fread(data.table = F, input = "./Data_Freezes/V1/snRNA/Tumor_Cell_Reclustered/Paths_TumorCellOnlyReclustered_SeuratObject.20201119.v1.tsv")
## input the barcode-manualsubcluster info
barcode2subclusterid_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/map_tumorsubclusterid/map_barcode_with_manual_tumorsubcluster_id/20201130.v1/Barcode2TumorSubclusterId.20201130.v1.tsv", data.table = F)
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20200902.v1/scrublet.run20200902_adj_cutoff.united_outputs.tsv", data.table = F)

# set parameters ----------------------------------------------------------
## set min.pct
logfc.threshold.run <- 0
min.pct.run <- 0.1
min.diff.pct.run <- 0.1
## spcify assay
assay_process <- "RNA"
cat(paste0("Assay: ", assay_process, "\n"))
cat("###########################################\n")


# for each aliquot, input seurat object and fetch data and write data --------------------
## merge cluster/scrublet info
barcode2info_df <- merge(x = barcode2subclusterid_df, y = barcode2scrublet_df, 
                         by.x = c("orig.ident", "barcode", "easy_id"), 
                         by.y = c("Aliquot", "Barcodes", "Aliquot_WU"), all.x = T)
markers_wilcox_df <- NULL
for (sample_tmp in unique(srat_paths_df$Sample)) {
  file2write <- paste0(dir_out,
                       sample_tmp, "TumorManualCluster.DEGs.Wilcox.Minpct", min.pct.run, 
                       ".Logfc", logfc.threshold.run,
                       ".min.diff.pct", min.diff.pct,
                       ".tsv")
  if (!file.exists(file2write)) {
    ## input seurat object
    seurat_obj_path <- srat_paths_df$Path_katmai[srat_paths_df$Sample == sample_tmp]
    seurat_obj_path
    srat <- readRDS(file = seurat_obj_path)
    ## change ident to manual subcluster
    ### get the barcode-subclusterid for this aliquot
    aliquot_tmp <- srat_paths_df$Aliquot.snRNA[srat_paths_df$Sample == sample_tmp]
    metadata_new_df <- barcode2info_df %>%
      filter(orig.ident == aliquot_tmp) %>%
      filter(!predicted_doublet)
    ## remove doublets
    srat <- subset(x = srat, cells = metadata_new_df$barcode)
    print(dim(srat))
    ### change meta data
    srat@meta.data$id_manual_cluster <- mapvalues(x = rownames(srat$meta.data), from = metadata_new_df$barcode, to = as.vector(metadata_new_df$id_manual_cluster_w0))
    ### change ident
    Idents(srat) <- "id_manual_cluster"
    
    ## find all markers using Wilcox testing
    markers_wilcox_tmp <- FindAllMarkers(object = srat, test.use = "wilcox", only.pos = F, 
                                         min.pct = min.pct.run, logfc.threshold = logfc.threshold.run, min.diff.pct = min.diff.pct.run, verbose = T, assay = assay_process)
    
    ## add the current DEGs into the super table
    if (nrow(markers_wilcox_tmp) > 0) {
      ## add cell numbers
      cellcounts_bycluster_df <- metadata_new_df %>%
        select(id_manual_cluster_w0) %>%
        table() %>%
        as.data.frame() %>%
        rename(id_manual_cluster = ".")
      cellcount_all <- sum(cellcounts_bycluster_df$Freq)
      cellcounts_bycluster_df <- cellcounts_bycluster_df %>%
        rename(cellcount_group1 = Freq) %>%
        mutate(cellcount_group2 = (cellcount_all - cellcount_group1))
      markers_wilcox_tmp <- merge(x = markers_wilcox_tmp, y = cellcounts_bycluster_df, by.x = c("cluster"), by.y = c("id_manual_cluster"), all.x = T)
      markers_wilcox_tmp$easy_id <- sample_tmp
      ## write output
      write.table(x = markers_wilcox_tmp, file = file2write, sep = "\t", quote = F, row.names = F)
      cat(paste0("finish writing for ", sample_tmp, "!\n\n"))
      
    }
  }
}


