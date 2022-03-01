# # library or install.packages-----------------------------------------------------------------
# packages = c(
#   "rstudioapi",
#   "D3GB",
#   "openxlsx",
#   "Matrix",
#   "optparse",
#   "bit64",
#   "Rtsne",
#   "Rmisc",
#   "ggplot2",
#   "RColorBrewer",
#   "ggrepel",
#   "ggthemes",
#   "data.table",
#   "pheatmap",
#   "cowplot",
#   'devtools', 
#   'readr',
#   'readxl',
#   'stringr',
#   'cowplot',
#   "Seurat",
#   'roxygen2'
# )
# 
# for (pkg_name_tmp in packages) {
#   if (!(pkg_name_tmp %in% installed.packages()[,1])) {
#     install.packages(pkg_name_tmp, dependencies = T)
#   }
# }
# 
# # install packages for plotting -------------------------------------------
# packages = c(
#   "rcartocolor",
#   "Polychrome",
#   "pals"
# )
# 
# for (pkg_name_tmp in packages) {
#   if (!(pkg_name_tmp %in% installed.packages()[,1])) {
#     install.packages(pkg_name_tmp, dependencies = T)
#   }
# }

# # install/library rjags -----------------------------------------------------------
# pkg_name_tmp <- "rjags"
# if (!(pkg_name_tmp %in% installed.packages()[,1])) {
#   devtools::install_url(url = "http://sourceforge.net/projects/mcmc-jags/files/rjags/4/rjags_4-4.tar.gz",
#                         args="--configure-args='--with-jags-include=/Users/casallas/homebrew/opt/jags/include/JAGS -with-jags-lib=/Users/casallas/homebrew/opt/jags/lib'")
# }
# library(package = pkg_name_tmp, character.only = T)
# 
# 
# 
# # install pkgs using BiocManager and library -----------------------------------------------------------
# packages = c(
#   "tibble",
#   "infercnv",
#   "ComplexHeatmap",
#   "circlize",
#   "R.methodsS3",
#   "irlba",
#   "MAST",
#   "SummarizedExperiment",
#   "DESeq2",
#   "rtracklayer",
#   "monocle",
#   "SingleCellExperiment",
#   "copynumber",
#   "dplyr"
# )
# 
# for (pkg_name_tmp in packages) {
#   if (!(pkg_name_tmp %in% installed.packages()[,1])) {
#     BiocManager::install(pkgs = pkg_name_tmp, update = F)
#   }
#   library(package = pkg_name_tmp, character.only = T)
# }


# ## install cellrangerRkit
# pkg_name_tmp <- "cellrangerRkit"
# if (!(pkg_name_tmp %in% installed.packages()[,1])) {
#   stop("run 'git clone https://github.com/hb-gitified/cellrangerRkit.git' on subdirectory called dependencies")
#   source(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/cellrangerRkit/scripts/rkit-install-2.0.0.R")
# }
# library(package = pkg_name_tmp, character.only = T)


# ## install loomR
# packages = c("loomR")
# 
# for (pkg_name_tmp in packages) {
#   if (!(pkg_name_tmp %in% installed.packages()[,1])) {
#     # install loomR from GitHub using the remotes package 
#     remotes::install_github(repo = 'mojaveazure/loomR', ref = 'develop')
#   }
#   library(package = pkg_name_tmp, character.only = T)
# }
# 
# 
# ## install Seurat
# pkg_name_tmp <- "Seurat"
# if (!(pkg_name_tmp %in% installed.packages()[,1])) {
#   install.packages(pkg_name_tmp, dependencies = T)
# }
# library(package = pkg_name_tmp, character.only = T)

# install SCENIC ----------------------------------------------------------
# ## Required
# BiocManager::install(c("AUCell", "RcisTarget"))
# BiocManager::install(c("GENIE3")) # Optional. Can be replaced by GRNBoost
# 
# ## Optional (but highly recommended):
# # To score the network on cells (i.e. run AUCell):
# BiocManager::install(c("zoo", "mixtools", "rbokeh"))
# # For various visualizations and perform t-SNEs:
# BiocManager::install(c("DT", "NMF", "pheatmap", "R2HTML", "Rtsne"))
# # To support paralell execution (not available in Windows):
# BiocManager::install(c("doMC", "doRNG"))

# if (!("SCENIC" %in% installed.packages()[,1])) {
#   packages = c("AUCell", "RcisTarget", "GENIE3", "zoo",  "mixtools", "rbokeh", "DT", "NMF", "pheatmap", "R2HTML", "Rtsne", "doMC", "doRNG")
#   for (pkg_name_tmp in packages) {
#     if (!(pkg_name_tmp %in% installed.packages()[,1])) {
#       BiocManager::install(pkgs = pkg_name_tmp, update = F)
#     }
#   }
#   if (!("SCopeLoomR" %in% installed.packages()[,1])) {
#     devtools::install_github("aertslab/SCopeLoomR", build_vignettes = TRUE)
#   }
#   
#   packageVersion("AUCell")
#   packageVersion("RcisTarget")
#   packageVersion("GENIE3")
#   
#   devtools::install_github("aertslab/SCENIC") 
#   packageVersion("SCENIC")
#   
# }
# 
# 
# # install TxDb.Hsapiens.UCSC.hg38.knownGene ----------------------------------------------------------
# if (!("TxDb.Hsapiens.UCSC.hg38.knownGene" %in% installed.packages()[,1])) {
#   BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
# }
# 
# 
# # install OmnipathR -------------------------------------------------------
# ## Last release in Bioconductor
# BiocManager::install("OmnipathR")
# ## Development version with the lastest updates
# BiocManager::install(version='devel')
# 
# 
# 
