# Yige Wu @WashU Sep 2019
## for integrating  snRNA datasets belong to the same patient

# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")

# set run id --------------------------------------------------------------
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)


# set output directory ----------------------------------------------------------
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# set case ids to be processed --------------------------------------------
case_id <- c("C3N-01200")

# input seurat processing summary ------------------------------------------------
seurat_summary <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/summary/ccRCC_snRNA_Downstream_Processing - Seurat_Preprocessing.20200116.v1.tsv", data.table = F)
seurat_summary2process <- seurat_summary %>%
  filter(Cellranger_reference_type == "pre-mRNA") %>%
  filter(Proceed_for_downstream == "Yes") %>%
  filter(Case %in% case_id) %>%
  filter(FACS == "") %>%
  mutate(Path_seurat_object = paste0("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/scRNA_auto/outputs/", Aliquot, FACS, 
                                     "/pf", `pre-filter.low.nCount`, "_fmin", low.nFeautre, "_fmax", high.nFeautre, "_cmin", low.nCount, "_cmax", high.nCount, "_mito_max", high.percent.mito, 
                                     "/", Aliquot, FACS, "_processed.rds"))
seurat_summary2process$Path_seurat_object

# set aliquot ids to be processed -----------------------------------------
snRNA_aliquot_ids <- seurat_summary2process$Aliquot

# create reference anchors  ----------------------------------------------------
renal.anchors <- readRDS(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/integration/integrate_same_patient_tumor_normal/20200117.v1/C3N-01200.Tummor_Segments.Anchors.20200117.v1.RDS")

# input marker table ------------------------------------------------------
gene2cellType_tab <- fread(input = "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Kidney_Markers/RCC_marker_gene_table_and_literature_review - Gene2CellType_Tab_w.HumanProteinAtlas.20200131.v1.tsv")

# Integrate Data -----------------------------------------------------------
## add cell type marker to the anchor gene set
new_anchor_features <- unique(c(renal.anchors@anchor.features, gene2cellType_tab$Gene))
length(new_anchor_features)
new_anchor_features %>% tail()

renal.integrated <- IntegrateData(anchorset = renal.anchors, dims = 1:20, features = new_anchor_features, features.to.integrate = new_anchor_features)
# Merging dataset 3 into 2
# Extracting anchors for merged samples
# Finding integration vectors
# Error in intI(i, n = d[1], dn[[1]], give.dn = FALSE) : 
#   invalid character indexing
rm(renal.anchors)

# switch to integrated assay. ---------------------------------------------
DefaultAssay(renal.integrated) <- "integrated" #only have 3000 features

# Run the standard workflow for visualization and clustering ------------
# The variable features of this assay are automatically set during IntegrateData
renal.integrated <- ScaleData(renal.integrated, verbose = F) 
renal.integrated <- RunPCA(renal.integrated, npcs = 30, verbose = FALSE)
renal.integrated <- RunUMAP(renal.integrated, reduction = "pca", dims = 1:20)
renal.integrated <- FindNeighbors(renal.integrated, reduction = "pca", dims = 1:20)
renal.integrated <- FindClusters(renal.integrated, resolution = 0.5)

saveRDS(object = renal.integrated, file = paste0(dir_out, case_id, ".Tummor_Normal.Integrated.", run_id, ".RDS"))
