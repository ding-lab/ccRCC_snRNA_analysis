# Yige Wu @WashU Aug 2020

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
path_rds <- "./Resources/Analysis_Results/recluster/recluster_cell_groups_in_integrated_data/subset_and_recluster/subset_stroma_and_recluster_on_katmai/20200824.v1/Stroma.Merged.20200824.v1.RDS"
srat <- readRDS(file = path_rds)
print("Finish reading RDS file\n")

# determine the statistical significance of PCA scores --------------------
# Do 200 random samplings to find significant genes, each time randomly permute 1% of genes
# This returns a 'p-value' for each gene in each PC, based on how likely the gene/PC score woud have been observed by chance
# Note that in this case we get the same result with 200 or 1000 samplings, so we do 200 here for expediency
srat <- JackStraw(object = srat, num.replicate = 100, dims = 20)
print("Finish running JackStraw\n")
## score the significance of the PCs
srat <- ScoreJackStraw(object = srat, dims = 1:20)
print("Finish running ScoreJackStraw\n")

# The jackStraw plot compares the distribution of P-values for each PC with a uniform distribution (dashed line)
# 'Significant' PCs will have a strong enrichment of genes with low p-values (solid curve above dashed line)
# In this case PC1-9 are strongly significant
file2write <- paste0(dir_out, "StromaReclustered.", "JackStrawPlot", ".png")
png(file2write, width = 1000, height = 1000, res = 150)
JackStrawPlot(pbmc, dims = 1:20)
print("Finish running JackStrawPlot\n")
dev.off()
print("Finish writing JackStrawPlot\n")

# identify genes significantly associated with any PC ---------------------
genes_sigpca = PCASigGenes(object = srat, pcs.use = 1:20, pval.cut = 1e-5, max.per.pc = 200)
print("Finish running PCASigGenes\n")
file2write <- paste0(dir_out, "StromaReclustered.", "PCASigGenes", run_id, ".RDS")
saveRDS(object = genes_sigpca, file = file2write, compress = T)
print("Finish writing PCASigGenes\n")



