# library or install.packages-----------------------------------------------------------------
packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  'stringr',
  "reshape2",
  "data.table",
  "Seurat",
  "Matrix",
  "optparse",
  "bit64",
  "Rtsne",
  "Rmisc",
  'devtools',
  'readr',
  'readxl',
  'cowplot',
  'roxygen2',
  "ggplot2",
  "RColorBrewer",
  "ggrepel",
  "ggthemes",
  "Polychrome",
  "ggrastr",
  "circlize"
)

for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    print(paste0(pkg_name_tmp, "is being installed!"))
    install.packages(pkg_name_tmp, dependencies = T)
  }
  print(paste0(pkg_name_tmp, " is installed!"))
}

# install pkgs using BiocManager and library -----------------------------------------------------------
packages = c(
  "ComplexHeatmap"
)

for (pkg_name_tmp in packages) {
  print(paste0(pkg_name_tmp, "is being installed!"))
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    BiocManager::install(pkgs = pkg_name_tmp, update = F)
  }
  print(paste0(pkg_name_tmp, " is installed!"))
}


# others ------------------------------------------------------------------
## install from https://www.xquartz.org/
