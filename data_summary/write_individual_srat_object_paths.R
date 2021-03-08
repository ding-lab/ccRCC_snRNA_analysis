# Yige Wu @WashU Jul 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input seurat processing info to input individual seurat object later
seurat_summary <- fread(input = "./Resources/snRNA_Processed_Data/scRNA_auto/summary/Seurat_Preprocessing.20210305.tsv", data.table = F)
seurat_summary2process <- seurat_summary %>%
  filter(Cellranger_reference_type == "pre-mRNA") %>%
  filter(Proceed_for_downstream == "Yes") %>%
  mutate(Path_box_seurat_object = paste0(dir_base, "Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                     "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                     "/", Aliquot, FACS, "_processed.rds")) %>%
  mutate(Path_katmai_seurat_object = gsub(x = Path_box_seurat_object, pattern = dir_base, replacement = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/")) %>%
  select(Case, Aliquot, Sample_Type, No.valid_barcode, Path_box_seurat_object, Path_katmai_seurat_object)

write.table(seurat_summary2process, file = paste0(dir_out, "Seurat_Object_Paths.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)
