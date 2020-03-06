# library or install.packages-----------------------------------------------------------------
packages = c(
  "rstudioapi",
  "openxlsx",
  "Matrix",
  "optparse",
  "bit64",
  "Rtsne",
  "Rmisc",
  "ggplot2",
  "RColorBrewer",
  "ggrepel",
  "data.table",
  "pheatmap",
  "cowplot",
  'devtools', 
  'readr',
  'readxl',
  'stringr',
  'cowplot',
  'roxygen2'
)

for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    install.packages(pkg_name_tmp, dependencies = T)
  }
  library(package = pkg_name_tmp, character.only = T)
}

# install/library rjags -----------------------------------------------------------
pkg_name_tmp <- "rjags"
if (!(pkg_name_tmp %in% installed.packages()[,1])) {
  devtools::install_url(url = "http://sourceforge.net/projects/mcmc-jags/files/rjags/4/rjags_4-4.tar.gz",
                        args="--configure-args='--with-jags-include=/Users/casallas/homebrew/opt/jags/include/JAGS -with-jags-lib=/Users/casallas/homebrew/opt/jags/lib'")
}
library(package = pkg_name_tmp, character.only = T)

# install pkgs using BiocManager and library -----------------------------------------------------------
packages = c(
  "tibble",
  "infercnv",
  "ComplexHeatmap",
  "circlize",
  "R.methodsS3",
  "irlba",
  "MAST",
  "SummarizedExperiment",
  "DESeq2",
  "rtracklayer",
  "monocle",
  "SingleCellExperiment",
  "dplyr"
)

for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    BiocManager::install(pkgs = pkg_name_tmp, update = F)
  }
  library(package = pkg_name_tmp, character.only = T)
}

# ## install cellrangerRkit
# pkg_name_tmp <- "cellrangerRkit"
# if (!(pkg_name_tmp %in% installed.packages()[,1])) {
#   stop("run 'git clone https://github.com/hb-gitified/cellrangerRkit.git' on subdirectory called dependencies")
#   source(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/cellrangerRkit/scripts/rkit-install-2.0.0.R")
# }
# library(package = pkg_name_tmp, character.only = T)


## install loomR
packages = c("loomR")

for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    # install loomR from GitHub using the remotes package 
    remotes::install_github(repo = 'mojaveazure/loomR', ref = 'develop')
  }
  library(package = pkg_name_tmp, character.only = T)
}


## install Seurat
pkg_name_tmp <- "Seurat"
if (!(pkg_name_tmp %in% installed.packages()[,1])) {
  install.packages(pkg_name_tmp, dependencies = T)
}
library(package = pkg_name_tmp, character.only = T)
