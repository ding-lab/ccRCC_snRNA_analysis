# Yige Wu @WashU Sep 2020
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
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir_katmai(path_this_script), run_id, "/")
dir.create(dir_out)
library(future)
plan("multiprocess", workers = 5)
options(future.globals.maxSize = 5 * 1024^3) # for 5 Gb RAM

# input dependencies ------------------------------------------------------
## input the integrated data
path_rds <- "./Resources/Analysis_Results/merging/33_aliquot_merging_without_anchoring/20210428.v2/33_aliquot_merged_without_anchoring.20210428.v2.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")
## input the barcode-cell-type table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups_33aliquots/20210423.v1/33Aliquot.Barcode2CellType.20210423.v1.tsv", data.table = F)
cat("finish reading the barcode-to-cell type table!\n")
## input idemta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210423.v1/meta_data.20210423.v1.tsv")
## input PBRM1 and BAP1 classification
mut_df <- fread(data.table = F, input = "./Resources/Analysis_Results/bulk/mutation/annotate_cptac_sample_by_pbrm1_bap1_mutation/20210310.v1/PBRM1_BAP1_Mutation_Status_By_Case.20210310.v1.tsv")

# set parameters for findmarkers ------------------------------------------
logfc.threshold.run <- 0
min.pct.run <- 0.1
min.diff.pct.run <- 0
## spcify assay
assay <- "RNA"
DefaultAssay(srat) <- assay
cat(paste0("Assay: ", assay, "\n"))
cat("###########################################\n")
## specify test
test_process <- "LR"
## specify cell groups to compare
ident.use.1 <- "BAP1-mutated Tumor cells"
ident.use.2 <- "PBRM1/Non-mutant Tumor cells"

# preprocess the Seurat object meta data---------------------------------------------
## get aliquot ids for the two groups
cases_group1 <- mut_df$Case[mut_df$mutation_category_sim %in% c("Both mutated", "BAP1 mutated")]; cases_group1 <- cases_group1[!(cases_group1 %in% c("C3L-01287"))]
cases_group1 <- cases_group1[cases_group1 %in% idmetadata_df$Case[idmetadata_df$snRNA_available]] ## 5 BAP1-mutated cases, 2 BAP1 & PBRM1 mutated cases
aliquots_group1 <- idmetadata_df$Aliquot.snRNA[idmetadata_df$snRNA_available & (idmetadata_df$Case %in% cases_group1) & idmetadata_df$Sample_Type == "Tumor"]
aliquots_group1 ## 10 samples
cases_group2 <- mut_df$Case[mut_df$mutation_category_sim %in% c("Non-mutants", "PBRM1 mutated")]; cases_group2 <- cases_group2[!(cases_group2 %in% c("C3L-00359"))]
aliquots_group2 <- idmetadata_df$Aliquot.snRNA[idmetadata_df$snRNA_available & idmetadata_df$Case %in% cases_group2 & idmetadata_df$Sample_Type == "Tumor"]
aliquots_group2 ## 19 samples
BC <- srat@meta.data %>% rownames
## get original barcode
srat@meta.data$original_barcode <- BC %>% strsplit("_") %>% lapply("[[",1) %>% unlist
head(srat@meta.data$original_barcode)
cat("finish adding the simple barcode!\n")
## make combined id for the seurat meta data
srat@meta.data$id_aliquot_barcode <- paste0(srat@meta.data$orig.ident, "_", srat@meta.data$original_barcode)
head(srat@meta.data$id_aliquot_barcode)
cat("finish adding unique id for each barcode in the seurat object!\n")
## make combined id for the barcode2celltype table
barcode2celltype_df <- barcode2celltype_df %>%
  mutate(id_aliquot_barcode = paste0(orig.ident, "_", individual_barcode)) %>%
  mutate(group_findmarkers = ifelse((Cell_group5 == "Tumor cells") & (orig.ident %in% aliquots_group1), 
                                    ident.use.1, 
                                    ifelse((Cell_group5 == "Tumor cells") & (orig.ident %in% aliquots_group2), 
                                           ident.use.2, 
                                           "Others")))
## map group label
srat@meta.data$group_findmarkers <- mapvalues(x = srat@meta.data$id_aliquot_barcode, from = barcode2celltype_df$id_aliquot_barcode, to = as.vector(barcode2celltype_df$group_findmarkers), warn_missing = F)
cat("finish adding group labels\n")
print(table(srat@meta.data$group_findmarkers))
Idents(srat) <- "group_findmarkers" 

# run findallmarkers by tumor------------------------------------------------------
## run findmarkers
deg_df <- FindMarkers(object = srat, test.use = test_process, ident.1 = ident.use.1, ident.2 = ident.use.2, only.pos = F,
                      min.pct = min.pct.run, logfc.threshold = logfc.threshold.run, min.diff.pct = min.diff.pct.run, verbose = T)
cat("finish FindMarkers\n")
deg_df$genesymbol_deg <- rownames(deg_df)
deg_df$cellnumber_tumorcells <- length(which(srat@meta.data$group_findmarkers == ident.use.1))
deg_df$cellnumber_ptcells <- length(which(srat@meta.data$group_findmarkers == ident.use.2))

## write output
file2write <- paste0(dir_out, test_process, 
                     ".logfc.threshold", logfc.threshold.run, 
                     ".min.pct", min.pct.run,
                     ".min.diff.pct", min.diff.pct.run,
                     ".Assay", assay,
                     ".tsv")
write.table(x = deg_df, file = file2write, sep = "\t", quote = F, row.names = F)
cat("finish writing the result!\n\n")



