# Yige Wu @WashU Mar 2022
## https://yoseflab.github.io/VISION/articles/web_only/Seurat.html
## before run script, run conda activate vision

# set up libraries and output directory -----------------------------------
dir_base = "~/Library/CloudStorage/Box-Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA"
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
path_rds <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/Analysis_Results/integration/seuratintegrate_34_ccRCC_samples/reciprocalPCA_integrate_34_ccRCC_samples/20220222.v1/ccRCC.34samples.SeuratIntegrated.20220222.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading the RDS file!\n")
## set paths to signature objects
signatures <- c("./Resources/Knowledge/Databases/MSigDB/msigdb_v7.4_GMTs/h.all.v7.4.symbols.gmt",
                "./Resources/Knowledge/Databases/MSigDB/msigdb_v7.4_GMTs/c2.cp.v7.4.symbols.gmt",
                "./Resources/Knowledge/Databases/MSigDB/msigdb_v7.4_GMTs/c5.go.bp.v7.4.symbols.gmt")

# process -----------------------------------------------------------------
# ## doing this otherwise will has the following message
# ### Dropping 'orig.ident' from meta data as it is of type 'character' and has more than 20 unique values.  If you want to include this meta data variable, convert it to a factor before providing the data frame to Vision
# srat@meta.data$orig.ident <- factor(srat@meta.data$orig.ident)
vision.obj <- Vision(srat, signatures = signatures, pool = F)
print("Finish creating the vision object!\n")
# Set the number of threads when running parallel computations
options(mc.cores = 15)
vision.obj <- analyze(vision.obj)
print("Finish analyze the vision object!\n")
file2write <- paste0(dir_out, "ccRCC.34samples.SeuratIntegrated.Vision.", run_id, ".RDS")
saveRDS(object = vision.obj, file = file2write, compress = T)

sigScores <- getSignatureScores(vision.obj)
file2write <- paste0(dir_out, "ccRCC.34samples.SeuratIntegrated.Vision.scores.", run_id, ".RDS")
saveRDS(object = sigScores, file = file2write, compress = T)
print("Finish getSignatureScores!\n")

sigCorr <- getSignatureAutocorrelation(vis_obj)
sigCorr$gene_set <- rownames(sigCorr)
print("Finish getSignatureAutocorrelation!\n")
file2write <- paste0(dir_out, "ccRCC.34samples.SeuratIntegrated.Vision.SignatureAutocorrelation.", run_id, ".tsv")
write.table(x = sigCorr, file = file2write, quote = F, sep = "\t", row.names = F)




