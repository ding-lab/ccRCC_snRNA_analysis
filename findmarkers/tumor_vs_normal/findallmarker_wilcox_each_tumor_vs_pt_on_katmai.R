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
dir_out <- "./Resources/snRNA_Processed_Data/Differentially_Expressed_Genes/Tumorcells_vs_PTcells/"

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

# set parameters for findmarkers ------------------------------------------
logfc.threshold.run <- 0
min.pct.run <- 0.1
min.diff.pct.run <- 0.1
## spcify assay
assay_process <- "RNA"
cat(paste0("Assay: ", assay_process, "\n"))
cat("###########################################\n")
## specify cell groups to compare
group1_findmarkers <- "Tumor cells"
group2_findmarkers <- "Proximal tubule cells from NATs"

# preprocess ----------------------------------------
aliquots_tumor <- idmetadata_df$Aliquot.snRNA[idmetadata_df$snRNA_available & idmetadata_df$Sample_Type == "Tumor"]
aliquots_tumor
aliquots_nat <- idmetadata_df$Aliquot.snRNA[idmetadata_df$snRNA_available & idmetadata_df$Sample_Type == "Normal"]
aliquots_nat

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
for (aliquot_tumor_tmp in aliquots_tumor) {
  easyid_tumor_tmp <- idmetadata_df$Aliquot.snRNA.WU[idmetadata_df$Aliquot.snRNA == aliquot_tumor_tmp]
  file2write <- paste0(dir_out,
                       easyid_tumor_tmp, ".Tumorcells_vs_PTcells.tsv")
  if (!file.exists(file2write)) {
    cat(paste0("Processing ", easyid_tumor_tmp, "\n"))
    ## make combined id for the barcode2celltype table
    barcode2celltype_df <- barcode2celltype_df %>%
      mutate(id_aliquot_barcode = paste0(orig.ident, "_", individual_barcode)) %>%
      mutate(group_findmarkers = ifelse(Cell_group5 == "Tumor cells" & orig.ident == aliquot_tumor_tmp, 
                                        "Tumor cells", 
                                        ifelse(Cell_type.detailed == "Proximal tubule" & orig.ident %in% aliquots_nat, "Proximal tubule cells from NATs", "Others")))
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
    
    ## write output
    write.table(x = deg_df, file = file2write, sep = "\t", quote = F, row.names = F)
    cat("finish writing the result!\n\n")
  }
}



