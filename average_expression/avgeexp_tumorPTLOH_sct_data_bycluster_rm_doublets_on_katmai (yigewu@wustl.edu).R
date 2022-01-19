# Yige Wu @WashU Jul 2020

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
## input the barcode-to-tumorsubcluster table
barcode2cellgroup_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_by_tumorcellmanualcluster_PTLOHclustes/20210903.v1/Barcode_byEpithelialCelltypes_BySample.20210903.v1.tsv", data.table = F)
# nrow(barcode2cellgroup_df)
cat("finish reading the barcode-to-tumorsubcluster table!\n")
barcode2scrublet_df <- fread(input = "./Resources/Analysis_Results/doublet/unite_scrublet_outputs/20210729.v1/scrublet.united_outputs.20210729.v1.tsv", data.table = F)
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210423.v1/meta_data.20210423.v1.tsv")
## input the integrated data (doublets removed)
path_rds <- "./Resources/Analysis_Results/merging/merge_35_samples/20210802.v1/RCC.35samples.Merged.20210802.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")
# spcify assay
assay_process <- "SCT"
slot_process <- "data"
cat(paste0("Assay: ", assay_process, "\n"))

# modify srat object ------------------------------------------------------
## add unique barcode id to the seurat meta data
### extract the meta data 
metadata_df <- srat@meta.data
metadata_df$barcode_31integration <- rownames(metadata_df)
metadata_df <- metadata_df %>%
  mutate(barcode_individual = str_split_fixed(string = barcode_31integration, pattern = "_", n = 2)[,1]) %>%
  mutate(aliquot_barcode = paste0(orig.ident, "_", barcode_individual))
### modify the tumor barcode table
barcode2scrublet_df <- barcode2scrublet_df %>%
  mutate(aliquot_barcode = paste0(Aliquot, "_", Barcode)) %>%
  filter(predicted_doublet)
### add new column for tumor cluster name
tmp <- mapvalues(x = metadata_df$aliquot_barcode, from = barcode2cellgroup_df$sample_barcode, to = as.vector(barcode2cellgroup_df$cell_group))
tmp[tmp %in% metadata_df$aliquot_barcode] <- "other"
tmp[tmp %in% barcode2scrublet_df$aliquot_barcode] <- "other"
table(tmp)
## change ident
srat@meta.data$Name_Subcluster <- tmp
Idents(srat) <- "Name_Subcluster"

## run average expression
aliquot.averages <- AverageExpression(srat, assays = assay_process, slot = slot_process)
print("Finish running AverageExpression!\n")
cat("###########################################\n")

## write output
file2write <- paste0(dir_out, "AverageExpression_ByTumorPTLOHSubcluster.", run_id, ".tsv")
write.table(aliquot.averages, file = file2write, quote = F, sep = "\t", row.names = T)
cat("Finished saving the output\n")
cat("###########################################\n")
