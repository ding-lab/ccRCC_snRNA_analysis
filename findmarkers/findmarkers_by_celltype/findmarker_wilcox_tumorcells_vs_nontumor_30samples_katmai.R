# Yige Wu @WashU Apr 2020
## finding differentially expressed gene for each cell type using integrared object

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
setwd(dir_base)
## load libraries
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "Seurat",
  "future",
  "future.apply"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
# set up future for parallization
plan("multiprocess", workers = 4)
options(future.globals.maxSize = 10000 * 1024^2)
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
source("./ccRCC_snRNA_analysis/functions.R")
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the integrated data
path_rds <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/Analysis_Results/merging/merge_34_ccRCC_samples/20220218.v1/ccRCC.34samples.Merged.20220218.v1.RDS"
srat <- readRDS(file = path_rds)
## input cell type per barcode table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups_35aliquots/20210802.v1/35Aliquot.Barcode2CellType.20210802.v1.tsv", data.table = F)
print("Finish reading the barcode2celltype_df file!\n")
## input meta data
metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210809.v1/meta_data.20210809.v1.tsv")

# set parameters for findmarkers ------------------------------------------
logfc.threshold.run <- 0.5
min.pct.run <- 0.1
min.diff.pct.run <- 0.1
DefaultAssay(srat)<-"RNA"
idents_group1 <- "Tumor cells"
idents_group2 <- c("Immune", "Normal epithelial cells", "Stroma")

# preprocess --------------------------------------------------------------
orig.idents_process <- metadata_df$Aliquot.snRNA[metadata_df$Sample_Type == "Tumor" & metadata_df$snRNA_available & metadata_df$Case != "C3L-00359"]

# set ident ---------------------------------------------------------------
## process barcode table
barcode2celltype_df <- barcode2celltype_df %>%
  mutate(cell_id = paste0(orig.ident, "_", individual_barcode))
## process seurat object
srat@meta.data$barcode_original <- str_split_fixed(string = rownames(srat@meta.data), pattern = "_", n = 2)[,1]
head(srat@meta.data$barcode_original)
srat@meta.data$cell_id <- paste0(srat@meta.data$orig.ident, "_", srat@meta.data$barcode_original)
srat@meta.data$Cell_group5 <- mapvalues(x = srat@meta.data$cell_id, from = barcode2celltype_df$cell_id, to = as.vector(barcode2celltype_df$Cell_group5))
srat@meta.data$Cell_group_test <- paste0(srat@meta.data$orig.ident, "_", srat@meta.data$Cell_group5)
Idents(srat) <- "Cell_group_test"
table(Idents(srat))
idents_group1 <- paste0(aliquots_group1, "_", idents_group1); idents_group1 <- idents_group1[idents_group1 %in% unique(Idents(srat))]; idents_group1
idents_group2 <- paste0(aliquots_group2, "_", idents_group2); idents_group2 <- idents_group2[idents_group2 %in% unique(Idents(srat))]; idents_group2

# run findallmarkers ------------------------------------------------------
markers <- FindMarkers(object = srat, test.use = "wilcox", ident.1 = idents_group1, ident.2 = idents_group2, only.pos = F,
                                min.pct = min.pct.run, logfc.threshold = logfc.threshold.run, min.diff.pct = min.diff.pct.run, verbose = T)
markers$gene_symbol <- rownames(markers)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "logfcthreshold.", logfc.threshold.run, ".minpct.", min.pct.run, ".mindiffpct.", min.diff.pct.run, ".tsv")
write.table(x = markers, file = file2write, quote = F, sep = "\t", row.names = F)


