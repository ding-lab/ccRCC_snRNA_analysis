# Yige Wu @WashU Sep 2020

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
## input the integrated data
path_rds <- "./Resources/snRNA_Processed_Data/scRNA_auto/outputs/CPT0001260013/pf1000_fmin200_fmax10000_cmin1000_cmax10000_mito_max0.1/CPT0001260013_processed.rds"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")
## input doublet prediction
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20210729.v1/scrublet.united_outputs.20210729.v1.tsv", data.table = F)
## input barcode-to-cell type info
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups_33aliquots/20210423.v1/33Aliquot.Barcode2CellType.20210423.v1.tsv")

# set parameters for findmarkers ------------------------------------------
easyid_process <- "C3L-00079-T1"
aliquot_process <- unique(barcode2scrublet_df$Aliquot[barcode2scrublet_df$Aliquot_WU == easyid_process])
barcode2celltype_df$Cell_group14_w_transitional[barcode2celltype_df$orig.ident == aliquot_process] %>% table()
## parameters for the findmarkers
# logfc.threshold.run <- log(2)
logfc.threshold.run <- 0
min.pct.run <- 0.1
min.diff.pct.run <- 0.1
## spcify assay
assay_process <- "RNA"
cat(paste0("Assay: ", assay_process, "\n"))
cat("###########################################\n")

# prepare seurat object ---------------------------------------------------
## subset
scrublets_df <- barcode2scrublet_df %>%
  filter(Aliquot_WU == easyid_process) %>%
  filter(predicted_doublet)
barcodes_keep <- rownames(srat@meta.data)
barcodes_keep <- barcodes_keep[!(barcodes_keep %in% scrublets_df$Barcodes)]
srat <- subset(x = srat, cells = barcodes_keep)
## add cell type
barcode2celltype_tmp_df <- barcode2celltype_df %>%
  filter(orig.ident == aliquot_process)
srat@meta.data$Cell_group14_w_transitional <- mapvalues(x = rownames(srat@meta.data), from = barcode2celltype_tmp_df$individual_barcode, to = as.vector(barcode2celltype_tmp_df$Cell_group14_w_transitional))
table(srat@meta.data$Cell_group14_w_transitional)
Idents(srat) <- "Cell_group14_w_transitional"

# run findallmarkers ------------------------------------------------------
degs_df <- FindMarkers(object = srat, test.use = "wilcox", ident.1 = "EMT tumor cells", ident.2 = "Tumor cells", only.pos = F,
                             min.pct = min.pct.run, logfc.threshold = logfc.threshold.run, min.diff.pct = min.diff.pct.run, verbose = T)
degs_df$row_name <- rownames(degs_df)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "findmarkers.wilcox.EMT_vs_nonEMT_tumorcells", 
                     ".logfcthreshold", logfc.threshold.run, ".minpct", min.pct.run, ".mindiffpct", min.diff.pct.run, ".tsv")
write.table(x = degs_df, file = file2write, sep = "\t", quote = F, row.names = F)
cat(paste0("Finished writing the output", "\n"))

