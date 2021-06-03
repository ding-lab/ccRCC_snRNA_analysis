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
plan("multiprocess", workers = 10)
options(future.globals.maxSize = 5 * 1024^3) # for 5 Gb RAM


# input dependencies ------------------------------------------------------
## input the integrated data
path_rds <- "./Data_Freezes/V2/snRNA/All_Cells_Merged/33_aliquot_merged_without_anchoring.20210428.v2.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file")
## input the barcode-cell-type table
barcode2celltype_df <- fread(input = "./Data_Freezes/V2/snRNA/Cell_Type_Assignment/33Aliquot.Barcode2CellType.20210423.v1.tsv", data.table = F)
cat("finish reading the barcode-to-cell type table!\n")
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20210423.v1/meta_data.20210423.v1.tsv")

PrctCellExpringGene <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(lapply(genes,calc_helper, object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = Seurat::SplitObject(object, group.by)
    factors = names(list)
    
    results = lapply(list, PrctCellExpringGene, genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}

calc_helper <- function(object,genes){
  counts = object[['RNA']]@counts
  ncells = ncol(counts)
  if(genes %in% row.names(counts)){
    sum(counts[genes,]>0)/ncells
  }else{return(NA)}
}

PrctCellExpringGeneParallel <- function(object, genes, group.by = "all"){
  if(group.by == "all"){
    prct = unlist(BiocParallel::bplapply(genes,calc_helper, BPPARAM = BiocParallel::MulticoreParam(log = T, workers = 4), object=object))
    result = data.frame(Markers = genes, Cell_proportion = prct)
    return(result)
  }
  
  else{        
    list = Seurat::SplitObject(object, group.by)
    factors = names(list)
    
    results = BiocParallel::bplapply(list, PrctCellExpringGene, BPPARAM = BiocParallel::MulticoreParam(log = T, workers = 4), genes=genes)
    for(i in 1:length(factors)){
      results[[i]]$Feature = factors[i]
    }
    combined = do.call("rbind", results)
    return(combined)
  }
}

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
## add easy id
srat@meta.data$id_easy <- mapvalues(x = srat@meta.data$orig.ident, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
barcode2celltype_df$id_easy <- mapvalues(x = barcode2celltype_df$orig.ident, from = idmetadata_df$Aliquot.snRNA, to = as.vector(idmetadata_df$Aliquot.snRNA.WU))
## make combined id for the barcode2celltype table
aliquots_nat <- idmetadata_df$Aliquot.snRNA[idmetadata_df$snRNA_available & idmetadata_df$Sample_Type == "Normal"]
aliquots_nat
barcode2celltype_df <- barcode2celltype_df %>%
  mutate(id_aliquot_barcode = paste0(orig.ident, "_", individual_barcode)) %>%
  mutate(group_findmarkers = ifelse(Cell_group5 == "Tumor cells", paste0("Tumorcells_", id_easy),
                                    ifelse(Cell_type.detailed == "Proximal tubule" & orig.ident %in% aliquots_nat, "Proximal tubule cells from NATs", "Others")))
## map group label
srat@meta.data$group_findmarkers <- mapvalues(x = srat@meta.data$id_aliquot_barcode, from = barcode2celltype_df$id_aliquot_barcode, to = as.vector(barcode2celltype_df$group_findmarkers), warn_missing = F)
cat("finish adding group labels\n")
print(table(srat@meta.data$group_findmarkers))
Idents(srat) <- "group_findmarkers" 

# process pct.expressed ---------------------------------------------------
pct_df <- PrctCellExpringGeneParallel(object = srat ,genes = rownames(srat@assays$RNA@counts), group.by = "group_findmarkers")
print("PrctCellExpringGene!")

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Pct.Expressed.byTumor_and_PT.", run_id, ".tsv")
write.table(x = pct_df, file = file2write, sep = '\t', quote = F, row.names = F)
print("Finished!")

