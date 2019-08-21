# Yige Wu @ WashU 2018 Jan
## shared parameters and functions for cptac2p_analysis


# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)
dir2dinglab_projects <- paste0(baseD, "Ding_Lab/Projects_Current/")
dir2cptac_pgdac <- paste0(dir2dinglab_projects, "CPTAC/PGDAC/")

# gene lists --------------------------------------------------------------
Plasma_markers <- data.frame(gene_symbol = c("SDC1", "IGHG1", "IGHG3", "IGHG4"), cell_type = "Plasma")
B_markers <- data.frame(gene_symbol = c("CD19", "MS4A1", "CD79A", "CD79B"), cell_type = "B")
Monocyte_markers <- data.frame(gene_symbol = c("LYZ", "CD14", "S100A8"), cell_type = "Monocyte")
Macrophage_markers <- data.frame(gene_symbol = c("FCGR3A", "MS4A7", "IFITM3"), cell_type = "Macrophage")
DC_markers <- data.frame(gene_symbol = c("FCER1A"), cell_type = "DC")
CD8T_markers <- data.frame(gene_symbol = c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B"), cell_type = "CD8+T")
CD4T_markers <- data.frame(gene_symbol = c("CD3D", "CD3E", "CD3G", "IL7R", "LDHB", "NOSIP", "CD4"), cell_type = "CD4+T")
NK_markers <- data.frame(gene_symbol = c("GNLY", "TYROBP", "HOPX", "FCGR3A"), cell_type = "NK")
Erythrocytes_markers <- data.frame(gene_symbol = c("HBD", "GYPA", "HBA1", "HBA2", "CA1"), cell_type = "Erythrocyte")

SMGs <- list()
SMGs[["CCRCC"]] <- c("VHL", "PBRM1", "SETD2", "KDM5C", "PTEN", "BAP1", "MTOR", "TP53")

# library -----------------------------------------------------------------
packages = c(
  "Matrix",
  "bit64",
  "Rtsne",
  "Rmisc",
  "ggplot2",
  "data.table",
  "pheatmap",
  'devtools', 
  'readr',
  'readxl',
  'stringr',
  'roxygen2'
)

for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    install.packages(pkg_name_tmp, dependencies = T)
  }
  library(package = pkg_name_tmp, character.only = T)
}

packages = c(
  "R.methodsS3",
  "irlba",
  "MAST",
  "SummarizedExperiment",
  "DESeq2",
  "rtracklayer",
  "monocle",
  "dplyr",
  "SingleCellExperiment"
)

for (pkg_name_tmp in packages) {
  if (!(pkg_name_tmp %in% installed.packages()[,1])) {
    BiocManager::install(pkgs = pkg_name_tmp)
  }
  library(package = pkg_name_tmp, character.only = T)
}

## install cellrangerRkit
pkg_name_tmp <- "cellrangerRkit"
if (!(pkg_name_tmp %in% installed.packages()[,1])) {
  stop("run 'git clone https://github.com/hb-gitified/cellrangerRkit.git' on subdirectory called dependencies")
  source(file = "./Ding_Lab/Projects_Current/RCC/ccRCC_single_cell/ccRCC_single_cell_analysis/scRNA/dependencies/cellrangerRkit/scripts/rkit-install-2.0.0.R")
}
library(package = pkg_name_tmp, character.only = T)


## install loomR
# install loomR from GitHub using the remotes package 
remotes::install_github(repo = 'mojaveazure/loomR', ref = 'develop')

## install Seurat
pkg_name_tmp <- "Seurat"
if (!(pkg_name_tmp %in% installed.packages()[,1])) {
  install.packages(pkg_name_tmp, dependencies = T)
}
library(package = pkg_name_tmp, character.only = T)


# Set paths ---------------------------------------------------------------
dir2dinglab_projects <- paste0(baseD, "Ding_Lab/Projects_Current/")
dir2current_project <- paste0(dir2dinglab_projects, "RCC/ccRCC_single_cell/")
dir2analysis_results <- paste0(dir2current_project, "analysis_results/")

###########################################
######## FUNCTIONS
###########################################

makeOutDir = function() {
  folders <- strsplit(x = rstudioapi::getSourceEditorContext()$path, split = "\\/")[[1]]
  folder_num <- which(folders == "ccRCC_single_cell_analysis") + 1
  dir2analysis_resultsnow <- paste(strsplit(paste(folders[folder_num:length(folders)], collapse = "/"), split = "\\.")[[1]][1], sep = "/")
  dir2analysis_resultsnow <- paste0(dir2analysis_results, dir2analysis_resultsnow, "/")
  dir.create(dir2analysis_resultsnow)
  dir2analysis_resultsnow_son <- dir2analysis_resultsnow
  dirs2make <- NULL
  while (!dir.exists(dir2analysis_resultsnow_son)) {
    tmp <- strsplit(dir2analysis_resultsnow_son, split = "\\/")[[1]]
    dir2analysis_resultsnow_parent <-paste(tmp[-length(tmp)], collapse = "/")
    dir.create(dir2analysis_resultsnow_parent)
    dir.create(dir2analysis_resultsnow_son)
    dir.create(dir2analysis_resultsnow)
    if (!dir.exists(dir2analysis_resultsnow_son)) {
      dirs2make[length(dirs2make) + 1] <- dir2analysis_resultsnow_son
    }
    dir2analysis_resultsnow_son <- dir2analysis_resultsnow_parent
  }
  
  if (length(dirs2make) > 0){
    for (i in 1:length(dirs2make)) {
      dir.create(dirs2make[i])
    }
  } 
  return(dir2analysis_resultsnow)
}

loadMaf <- function() {
  maf <- fread(input = paste0(dir2cptac_pgdac, "ccRCC_discovery_manuscript/ccRCC_expression_matrices/Somatic_Variants/ccrcc.somatic.consensus.gdc.umichigan.wu.112918.maf"), data.table = F, fill=TRUE) 
  maf <- data.frame(maf)
  print(paste0("MAF has ", nrow(maf), " lines\n"))
  return(maf)
}


loadCNAstatus <- function() {
  ## input CNA values
  cancer <- "CCRCC"
  cna <- fread(input = paste0("./Ding_Lab/Projects_Current/PanCan_Phospho-signaling/analysis_results/preprocess_files/tables/parse_", cancer, "_data_freeze/somatic_CNA.", cancer, ".partID.txt"), data.table = F)
  cna_head <- cna$gene
  cna_mat <- cna[, colnames(cna)[!(colnames(cna) %in% "gene")]]
  cna_status <- matrix(data = "neutral", nrow = nrow(cna_mat), ncol = ncol(cna_mat))
  cna_status[cna_mat > 0.1] <- "amplification"
  cna_status[cna_mat < -0.1] <- "deletion"
  cna_status <- data.frame(cbind(cna$gene, cna_status))
  colnames(cna_status) <- colnames(cna)
  return(cna_status)
}

loadRNA <- function() {
  rna_fn <- "ccRcc_RNA_rpkm_Mich_formatted_tumor.csv"
  rna_tab <- fread(input = paste0(dir2dinglab_projects, "TP53_shared_data/resources/rna/", rna_fn), data.table = F)
  
  colnames(rna_tab)[1] <- "gene"
  rna_mat <- as.matrix(rna_tab[,-1])
  rna_mat_log2 <- log2(rna_mat+1)
  rna_tab <- data.frame(gene = rna_tab$gene)
  rna_tab <- cbind(rna_tab, as.data.frame(rna_mat_log2))
  return(rna_tab) 
}

loadParseProteomicsData <- function(expression_type, sample_type) {
  ## expresson_type: PRO or PHO (phosphosite level) or PHO_collapsed (protein level)
  ## sample_type: tumor or normal
  ## pipeline_type: CDAP or PGDAC
  ## norm_type: unnormalized or scaled
  cancer <- "CCRCC"
  pipeline_type <-  "PGDAC"
  norm_type <- "MD_MAD"
  
  dir1 <- paste0("./Ding_Lab/Projects_Current/PanCan_Phospho-signaling/analysis_results/preprocess_files/tables/parse_", cancer, "_data_freeze", "" , "/")
  exp_data <- fread(input = paste0(dir1, cancer, "_", expression_type, "_", sample_type, "_", pipeline_type, "_", norm_type, "_", "partID", ".txt"), data.table = F)
  return(exp_data)
}

generate_somatic_mutation_matrix <- function(pair_tab, maf) {
  genes4mat <- unique(unlist(pair_tab))
  length(genes4mat)
  
  maf <- maf[maf$Hugo_Symbol %in% genes4mat,]
  nrow(maf)
  maf$sampID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
  
  mut_mat <- reshape2::dcast(data = maf, Hugo_Symbol ~ sampID, fun =  function(x) {
    variant_class <- paste0(unique(x), collapse = ",")
    return(variant_class)
  }, value.var = "Variant_Classification", drop=FALSE)
  rownames(mut_mat) <- as.vector(mut_mat$Hugo_Symbol)
  return(mut_mat)
}

scale_by_row = function(m){
  m = as.matrix(m)
  m[!is.finite(m)] = NA
  for (i in 1:nrow(m)){
    mean_tmp <-  mean(m[i,], na.rm=T)
    sd_tmp <- sd(m[i,], na.rm=T)
    if ( !is.nan(mean_tmp) & !is.na(sd_tmp) ) {
      m[i,] = (m[i,] - mean_tmp)/sd_tmp
    } else {
      m[i,] = NA
    }
  }
  return(m)
}

save_pheatmap_pdf <- function(x, filename, width=6, height=6) {
  pdf(filename, width = width, height = height)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png <- function(x, filename, width=1200, height=1000, res = 150) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}
