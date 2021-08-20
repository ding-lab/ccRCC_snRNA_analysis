# Yige Wu @WashU Jun 2021
## use RNA assay according to https://github.com/satijalab/seurat/issues/2646

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
options(future.globals.maxSize = 5 * 1024^3) # for 5 Gb RAM

# input dependencies ------------------------------------------------------
## input the integrated data
path_rds <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/Analysis_Results/merging/merge_35_samples/20210802.v1/RCC.35samples.Merged.20210802.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")
## input the barcode-cell-type table
barcode2celltype_df <- fread(input = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups_35aliquots/20210802.v1/35Aliquot.Barcode2CellType.20210802.v1.tsv", data.table = F)
cat("finish reading the barcode-to-cell type table!\n")
## input id data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210809.v1/meta_data.20210809.v1.tsv")
## input clinical data
clinical_bycase_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/clinical/extract_case_clinical_data/20201125.v1/snRNA_ccRCC_Clinicl_Table.20201125.v1.tsv")

# set parameters for findmarkers ------------------------------------------
logfc.threshold.run <- 0.25
min.pct.run <- 0.1
min.diff.pct.run <- 0.1
## spcify assay
assay_process <- "RNA"
DefaultAssay(srat) <- assay_process
cat(paste0("Assay: ", assay_process, "\n"))
cat("###########################################\n")
## specify test
test_process <- "wilcox"
## specify cell groups to compare
group1_findmarkers <- "Male Tumor cells"
group2_findmarkers <- "Female Tumor cells"
## specify case feature
case_feature_group1 <- "Male"
case_feature_group2 <- "Female"

# preprocess ----------------------------------------
cases_process_group1 <- clinical_bycase_df$Case[clinical_bycase_df$Gender == case_feature_group1]
aliquots_group1<- idmetadata_df$Aliquot.snRNA[idmetadata_df$snRNA_available & idmetadata_df$Sample_Type == "Tumor" & idmetadata_df$Case != "C3L-00359" & idmetadata_df$Case %in% cases_process_group1]
aliquots_group1
cases_process_group2 <- clinical_bycase_df$Case[clinical_bycase_df$Gender == case_feature_group2]
aliquots_group2<- idmetadata_df$Aliquot.snRNA[idmetadata_df$snRNA_available & idmetadata_df$Sample_Type == "Tumor" & idmetadata_df$Case != "C3L-00359" & idmetadata_df$Case %in% cases_process_group2]
aliquots_group2

# preprocess the Seurat object meta data---------------------------------------------
BC <- srat@meta.data %>% rownames
## get original barcode
srat@meta.data$original_barcode <- BC %>% strsplit("_") %>% lapply("[[",1) %>% unlist
head(srat@meta.data$original_barcode)
cat("finish adding the simple barcode!\n")
## make combined id for the seurat meta data
srat@meta.data$id_aliquot_barcode <- paste0(srat@meta.data$orig.ident, "_", srat@meta.data$original_barcode)
head(srat@meta.data$id_aliquot_barcode)
cat("finish adding unique id for each barcode in the seurat object!\n")

# run findallmarkers by tumor------------------------------------------------------
## make combined id for the barcode2celltype table
barcode2celltype_df <- barcode2celltype_df %>%
  mutate(id_aliquot_barcode = paste0(orig.ident, "_", individual_barcode)) %>%
  mutate(group_findmarkers = ifelse(Cell_group5 == "Tumor cells" & orig.ident %in% aliquots_group1, 
                                    group1_findmarkers, 
                                    ifelse(Cell_group5 == "Tumor cells" & orig.ident %in% aliquots_group2, group2_findmarkers, "Others")))
## map group label
srat@meta.data$group_findmarkers <- mapvalues(x = srat@meta.data$id_aliquot_barcode, from = barcode2celltype_df$id_aliquot_barcode, to = as.vector(barcode2celltype_df$group_findmarkers), warn_missing = F)
cat("finish adding group labels\n")

print(table(srat@meta.data$group_findmarkers))
Idents(srat) <- "group_findmarkers" 
## run findmarkers
deg_df <- FindMarkers(object = srat, test.use = test_process, ident.1 = group1_findmarkers, ident.2 = group2_findmarkers, only.pos = F,
                      min.pct = min.pct.run, logfc.threshold = logfc.threshold.run, min.diff.pct = min.diff.pct.run, verbose = T)
cat("finish FindMarkers\n")
deg_df$genesymbol_deg <- rownames(deg_df)
deg_df$cellnumber <- length(which(srat@meta.data$group_findmarkers == group1_findmarkers))
deg_df$cellnumber <- length(which(srat@meta.data$group_findmarkers == group2_findmarkers))

## write output
file2write <- paste0(dir_out, "Male_vs_Female", ".ccRCC.", test_process, 
                     ".logfc.threshold", logfc.threshold.run, 
                     ".min.pct", min.pct.run,
                     ".min.diff.pct", min.diff.pct.run,
                     ".Assay", assay_process,
                     ".tsv")
write.table(x = deg_df, file = file2write, sep = "\t", quote = F, row.names = F)
cat("finish writing the result!\n\n")





