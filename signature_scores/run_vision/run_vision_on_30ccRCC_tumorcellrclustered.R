# Yige Wu @WashU Mar 2022
## https://yoseflab.github.io/VISION/articles/web_only/Seurat.html
## before run script, run conda activate vision

# set up libraries and output directory -----------------------------------
dir_base = "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/"
setwd(dir_base)
packages = c(
  "Seurat",
  "VISION"
)
if (!("VISION" %in% installed.packages()[,1])) {
  print(paste0(pkg_name_tmp, "is being installed!"))
  devtools::install_github("YosefLab/VISION", dependencies = T)
}
for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
dir_out <- paste0("/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/Analysis_Results/signature_scores/run_vision/run_vision_on_34ccRCC_integrated/", run_id, "/")
dir.create(dir_out)

# input  ------------------------------------------------------
## input seurat object
path_rds <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/Analysis_Results/integration/seuratintegrate_34_ccRCC_samples/reciprocalPCA_integrate_30_ccRCC_tumorcells/20220404.v1/ccRCC.34samples.Tumorcells.SeuratIntegrated.20220404.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading the RDS file!\n")
## set paths to signature objects
signatures <- c("./Resources/Knowledge/Databases/MSigDB/msigdb_v7.4_GMTs/h.all.v7.4.symbols.gmt",
                "./Resources/Knowledge/Databases/MSigDB/msigdb_v7.4_GMTs/c2.cp.v7.4.symbols.gmt",
                "./Resources/Knowledge/Databases/MSigDB/msigdb_v7.4_GMTs/c5.go.bp.v7.4.symbols.gmt")

# process -----------------------------------------------------------------
## this is important because the default assay for the seurat object is "integrated" and I think if not set to RNA assay Vision will just take integrated assay,
## which will give an error about  "rownames(data) = NULL. Expression matrix must have gene names as the rownames" because
### > rownames(srat$integrated@counts)
### NULL
### rownames(srat$integrated@data will give top variably expressed genes
DefaultAssay(srat) <- "RNA"
vision.obj <- Vision(srat, signatures = signatures, pool = F)
print("Finish creating the vision object!\n")
# Set the number of threads when running parallel computations
options(mc.cores = 8)
vision.obj <- analyze(vision.obj)
print("Finish analyze the vision object!\n")
sigScores <- getSignatureScores(vision.obj)
print("Finish getSignatureScores!\n")

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "ccRCC.30ccRCC.TumorCellsReclustered.Vision.", run_id, ".RDS")
saveRDS(object = vision.obj, file = file2write, compress = T)
file2write <- paste0(dir_out, "ccRCC.30ccRCC.TumorCellsReclustered.Vision.sigScores.", run_id, ".RDS")
saveRDS(object = sigScores, file = file2write, compress = T)
