# Yige Wu @WashU Mar 2021
## for merging 32 snRNA datasets

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
# dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
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
options(future.globals.maxSize = 1000 * 1024^2)
library(harmony)

# input dependencies ------------------------------------------------------
path_rds <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/Analysis_Results/merging/merge_34_ccRCC_samples/20220218.v1/ccRCC.34samples.Merged.20220218.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading the RDS file!\n")
## input cell type per barcode table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups_35aliquots/20210802.v1/35Aliquot.Barcode2CellType.20210802.v1.tsv", data.table = F)
print("Finish reading the barcode2celltype_df file!\n")

# run harmony -------------------------------------------------------------
srat <- RunHarmony(object = srat, group.by.vars = "orig.ident", assay.use = "SCT")
cat("Finished RunHarmony!\n")
srat <- srat %>% 
  RunUMAP(reduction = "harmony", dims = 1:30, verbose = F) %>% 
  FindNeighbors(reduction = "harmony", k.param = 10, dims = 1:30) %>% 
  FindClusters() %>% 
  identity()
cat("Finished RunUMAP!\n")

## save as RDS file
file2write <- paste0(dir_out, "ccRCC.34samples.Merged.HarmonyIntegrated.20220218.v1.RDS")
saveRDS(object = srat, file = file2write, compress = T)
cat("Finished saving the output!\n")

# fetch data --------------------------------------------------------------
umap_data <- Seurat::FetchData(object = srat, vars = c("orig.ident", "ident", "UMAP_1", "UMAP_2"))
umap_data$barcode <- rownames(umap_data)

file2write <- paste0(dir_out, "ccRCC.34Sample.Merged.HarmonyIntegrated.Metadata.", run_id, ".tsv")
write.table(x = umap_data, file = file2write, quote = F, sep = "\t", row.names = F)
print("Finish writing the output!\n")

# plot --------------------------------------------------------------------
file2write <- paste0(dir_out, "ccRCC.34samples.Merged.HarmonyIntegrated.", run_id, ".png")
png(filename = file2write, width = 1000, height = 1100, res = 150)
DimPlot(srat,reduction = "umap")
dev.off()



