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
path_rds <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/Analysis_Results/merging/merge_34_ccRCC_samples/20220218.v1/ccRCC.34samples.Merged.20220218.v1.RDS"
srat <- readRDS(file = path_rds)
## input cell type per barcode table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups_35aliquots/20210802.v1/35Aliquot.Barcode2CellType.20210802.v1.tsv", data.table = F)
print("Finish reading the barcode2celltype_df file!\n")
## input meta data
metadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210809.v1/meta_data.20210809.v1.tsv")

# set parameters for findmarkers ------------------------------------------
logfc.threshold.run <- 0
min.pct.run <- 0.1
min.diff.pct.run <- 0
DefaultAssay(srat)<-"RNA"
celltypes_process <- unique(barcode2celltype_df$Cell_type.detailed[grepl(pattern = "Macrophage", x = barcode2celltype_df$Cell_type.detailed)])
celltypes_process <- c(celltypes_process, "TRM VSIR+")
samples_inflamm <- c("C3L-00416-T1", "C3L-01302-T1", "C3N-00437-T1", "C3N-01200-T2", "C3N-01200-T3")
samples_noninflamm <- c("C3L-00610-T1", "C3L-00813-T1", "C3L-00908-T1", "C3L-00917-T1", "C3L-01287-T1")

# preprocess --------------------------------------------------------------
barcode2celltype_df$sample_id <- mapvalues(x = barcode2celltype_df$orig.ident, from = metadata_df$Aliquot.snRNA, to = as.vector(metadata_df$Aliquot.snRNA.WU))
barcode2celltype_df <- barcode2celltype_df %>%
  mutate(id_aliquot_barcode = paste0(orig.ident, "_", individual_barcode)) %>%
  mutate(group_findmarkers = ifelse(sample_id %in% samples_inflamm & Cell_type.detailed %in% celltypes_process, "group1",
                                    ifelse(sample_id %in% samples_noninflamm& Cell_type.detailed %in% celltypes_process, "group2", "others")))
table(barcode2celltype_df$group_findmarkers)
print("Finished making cell group!")

# set ident ---------------------------------------------------------------
## make combined id for the Seurat meta data
srat@meta.data$id_aliquot_barcode <- paste0(srat@meta.data$orig.ident, "_", srat@meta.data$original_barcode)
srat@meta.data$group_findmarkers <- mapvalues(x = srat@meta.data$id_aliquot_barcode, from = barcode2celltype_df$id_aliquot_barcode, to = as.vector(barcode2celltype_df$group_findmarkers))
Idents(srat) <- "group_findmarkers"
table(Idents(srat))
# srat <- subset(srat, idents = c("group1", "group2"))
# table(Idents(srat))
# print("Finished subsetting!")

# run findallmarkers ------------------------------------------------------
deg_df <- FindMarkers(object = srat, test.use = "wilcox", ident.1 = "group1", ident.2 = "group2", only.pos = F,
                      min.pct = min.pct.run, logfc.threshold = logfc.threshold.run, min.diff.pct = min.diff.pct.run, verbose = T)
deg_df$gene_symbol <- rownames(deg_df)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "logfcthreshold.", logfc.threshold.run, ".minpct.", min.pct.run, ".mindiffpct.", min.diff.pct.run, ".tsv")
write.table(x = deg_df, file = file2write, quote = F, sep = "\t", row.names = F)


