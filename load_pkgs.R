# load libraries ----------------------------------------------------------

packages = c(
  "rstudioapi",
  "plyr",
  "dplyr",
  "stringr",
  "reshape2",
  "data.table",
  "Seurat"
)

for (pkg_name_tmp in packages) {
  library(package = pkg_name_tmp, character.only = T)
}