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
dir_parent_out <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/Analysis_Results/signature_scores/run_vision/getSignatureAutocorrelation_30ccRCC_tumorcellreclustered/"
dir.create(dir_parent_out)
dir_out <- paste0(dir_parent_out, run_id, "/")
dir.create(dir_out)

# input  ------------------------------------------------------
## input seurat object
path_rds <- "/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/Analysis_Results/signature_scores/run_vision/run_vision_on_30ccRCC_tumorcellreclustered/20220406.v1/ccRCC.30ccRCC.TumorCellsReclustered.Vision.20220406.v1.RDS"
vis_obj <- readRDS(file = path_rds)
print("Finish reading the RDS file!\n")

# process -----------------------------------------------------------------
options(mc.cores = 8)
sigCorr <- getSignatureAutocorrelation(vis_obj)
sigCorr$gene_set <- rownames(sigCorr)
print(head(sigCorr))
print("Finish getSignatureAutocorrelation!\n")

# save output -------------------------------------------------------------
file2write <- paste0(dir_out, "ccRCC.30ccRCC.TumorCellsReclustered.Vision.SignatureAutocorrelation.", run_id, ".tsv")
write.table(x = sigCorr, file = file2write, quote = F, sep = "\t", row.names = F)
