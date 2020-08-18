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
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the integrated data
path_rds <- "./Resources/Analysis_Results/integration/30_aliquot_integration/subset/subset_tumor_cells/20200309.v1/30_aliquot_integration.20200212.v3.RDS.tumor_cells.20200309.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")
cat("finish reading the barcode-to-tumorsubcluster table!\n")
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200505.v1/meta_data.20200505.v1.tsv")
## set aliquot ids for groups
aliquotids_group1 <- c("CPT0075720013", "CPT0063630004")
aliquotids_group2 <- unique(srat@meta.data$orig.ident)[!(unique(srat@meta.data$orig.ident) %in% c(aliquotids_group1))]

# set parameters for findmarkers ------------------------------------------
logfc.threshold.run <- 0.1
min.pct.run <- 0.1
min.diff.pct.run <- 0.1

# findmarkers -------------------------------------------------------------
markers_all_df <- NULL
for (aliquotid_group1 in aliquotids_group1) {
  for (aliquotid_group2 in aliquotids_group2) {
    cat(paste0("Starting ", aliquotid_group1, " vs ", aliquotid_group2, "###########################################\n"))
    # set ident ---------------------------------------------------------------
    srat@meta.data$group_findmarkers <- ifelse(srat@meta.data$orig.ident %in% aliquotid_group1, "group1", 
                                               ifelse(srat@meta.data$orig.ident %in% aliquotid_group2, "group2", "other"))
    Idents(srat) <- "group_findmarkers"

    markers_df <- FindMarkers(object = srat, ident.1 = "group1", ident.2 = "group2", test.use = "wilcox", only.pos = F, 
                              min.pct = min.pct.run, logfc.threshold = logfc.threshold.run, min.diff.pct = min.diff.pct.run, verbose = T)
    cat(paste0("Finished ", aliquotid_group1, " vs ", aliquotid_group2, "###########################################\n"))
    markers_df$gene <- rownames(markers_df)
    markers_df$ident.1 <- idmetadata_df$Aliquot.snRNA.WU[idmetadata_df$Aliquot.snRNA == aliquotid_group1]
    markers_df$ident.2 <- idmetadata_df$Aliquot.snRNA.WU[idmetadata_df$Aliquot.snRNA == aliquotid_group2]
    markers_all_df <- rbind(markers_all_df, markers_df)
  }
}
cat("###########################################\n")

# write output -------------------------------------------------------------
file2write <- paste0(dir_out, "findallmarkers_wilcox_1213_1302_vs_others.", "bysample.", ".logfcthreshold", logfc.threshold.run, ".minpct", min.pct.run, ".mindiffpct", min.diff.pct.run, ".tsv")
write.table(x = markers_all_df, file = file2write, sep = "\t", quote = F, row.names = F)





