# Yige Wu @WashU Sep 2020
## use RNA assay according to https://github.com/satijalab/seurat/issues/2646
## 2020-09-29 changed the sample to 7 snATAC tumors and 2 snATAC NATs with cell type 20200917.v2

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
path_rds <- "./Resources/Analysis_Results/integration/31_aliquot_integration/31_aliquot_integration_without_anchoring/20200727.v1/31_aliquot_integration_without_anchoring.20200727.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")
## input the barcode-cell-type table
barcode2celltype_df <- fread(input = "./Resources/Analysis_Results/annotate_barcode/annotate_barcode_with_major_cellgroups/20201121.v1/31Aliquot.Barcode2CellType.20201121.v1.tsv", data.table = F)
cat("finish reading the barcode-to-cell type table!\n")
## input idemta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200716.v1/meta_data.20200716.v1.tsv")

# set parameters for findmarkers ------------------------------------------
logfc.threshold.run <- 0
min.pct.run <- 0.1
min.diff.pct.run <- 0.1
## spcify assay
assay_process <- "RNA"
cat(paste0("Assay: ", assay_process, "\n"))
cat("###########################################\n")
## specify cell groups to compare
group1_findmarkers <- "Macrophages_Tumor"
group2_findmarkers <- "Macrophages_NAT"

# subset to cases with snATAC data ----------------------------------------
## Based on the small cell numbers of the C3L-00088-N and C3L-00088-T1, we exclude them
easyids_snatac <- c("C3N-00733-T1", "C3L-00610-T1", "C3L-01313-T1", "C3L-00416-T2", "C3L-01287-T1", "C3L-00917-T1", "C3N-01200-T1", "C3L-00088-T2", "C3L-00079-T1", "C3L-00448-T1",
                    "C3N-01200-N")
aliquots_snatac <- idmetadata_df$Aliquot.snRNA[idmetadata_df$Aliquot.snRNA.WU %in% easyids_snatac]
aliquots_snatac
aliquots_snatac_nat <- idmetadata_df$Aliquot.snRNA[idmetadata_df$Aliquot.snRNA.WU %in% easyids_snatac & idmetadata_df$Sample_Type == "Normal"]
aliquots_snatac_nat
aliquots_snatac_tumor <- idmetadata_df$Aliquot.snRNA[idmetadata_df$Aliquot.snRNA.WU %in% easyids_snatac & idmetadata_df$Sample_Type == "Tumor"]
aliquots_snatac_tumor
Idents(srat) <- "orig.ident"
srat <- subset(srat, idents = aliquots_snatac)
cat("finish subsetting the object!\n")

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
degs_combined_df <- NULL
for (aliquot_tumor_tmp in aliquots_snatac_tumor) {
  easyid_tumor_tmp <- idmetadata_df$Aliquot.snRNA.WU[idmetadata_df$Aliquot.snRNA == aliquot_tumor_tmp]
  cat(paste0("Processing ", easyid_tumor_tmp, "\n"))
  ## make combined id for the barcode2celltype table
  barcode2celltype_df <- barcode2celltype_df %>%
    mutate(id_aliquot_barcode = paste0(orig.ident, "_", individual_barcode)) %>%
    mutate(group_findmarkers = ifelse(Cell_group13 == "Macrophages" & orig.ident == aliquot_tumor_tmp, 
                                      "Macrophages_Tumor", 
                                      ifelse(Cell_type.detailed == "Macrophages" & orig.ident %in% aliquots_snatac_nat, "Macrophages_NAT", "Others")))
  
  ## map group label
  srat@meta.data$group_findmarkers <- mapvalues(x = srat@meta.data$id_aliquot_barcode, from = barcode2celltype_df$id_aliquot_barcode, to = as.vector(barcode2celltype_df$group_findmarkers), warn_missing = F)
  cat("finish adding group labels\n")
  
  print(table(srat@meta.data$group_findmarkers))
  Idents(srat) <- "group_findmarkers" 
  ## run findmarkers
  deg_df <- FindMarkers(object = srat, test.use = "wilcox", ident.1 = group1_findmarkers, ident.2 = group2_findmarkers, only.pos = F,
                        min.pct = min.pct.run, logfc.threshold = logfc.threshold.run, min.diff.pct = min.diff.pct.run, verbose = T)
  deg_df$genesymbol_deg <- rownames(deg_df)
  deg_df$easyid_tumor <- easyid_tumor_tmp
  deg_df$aliquot_tumor <- aliquot_tumor_tmp
  deg_df$cellnumber_tumorcells <- length(which(srat@meta.data$group_findmarkers == group1_findmarkers))
  deg_df$cellnumber_ptcells <- length(which(srat@meta.data$group_findmarkers == group2_findmarkers))
  
  ## combine with the super tale
  degs_combined_df <- rbind(degs_combined_df, deg_df)
}

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "findallmarkers_wilcox_each_snatac_tumor_vs_nat_macrophages.", run_id, ".tsv")
write.table(x = degs_combined_df, file = file2write, sep = "\t", quote = F, row.names = F)
cat("finish writing the result!\n")

