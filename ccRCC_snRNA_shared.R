# Yige Wu @ WashU 2018 Jan
## shared parameters and functions for cptac2p_analysis


# set working directory ---------------------------------------------------
baseD = "~/Box/"
setwd(baseD)

# Set paths ---------------------------------------------------------------
dir2dinglab_projects <- paste0(baseD, "Ding_Lab/Projects_Current/")
dir2current_project <- paste0(dir2dinglab_projects, "RCC/ccRCC_snRNA/")
dir2analysis_results <- paste0(dir2current_project, "Resources/Analysis_Results/")
dir2cptac_pgdac <- paste0(dir2dinglab_projects, "CPTAC/PGDAC/")

# load libraries ----------------------------------------------------------
packages = c(
  "rstudioapi",
  "optparse",
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

# set shared snRNA processing parameters ----------------------------------
num_pc <- 30


# gene lists --------------------------------------------------------------
## significantly mutated genes
ccRCC_SMGs <- c("VHL", 
                "PBRM1", 
                "SETD2", "BAP1", "KDM5C",   
                "PTEN", 
                "MTOR", "TSC1",
                "TP53")
## PBAF gens
### reference: https://www.nature.com/articles/onc20094/figures/1
### reference: https://www.nature.com/articles/onc20094/tables/1
pbaf_genes <- c("PBRM1", "ARID2", "SMARCE1", "SMARCC2", "ACTL6A", "ACTL6B", "SMARCC1", "SMARCD1", "SMARCB1")

# make output directory ---------------------------------------------------
makeOutDir = function() {
  folders <- strsplit(x = rstudioapi::getSourceEditorContext()$path, split = "\\/")[[1]]
  folder_num <- which(folders == "ccRCC_snRNA_analysis") + 1
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


# Copy number related functions and varaibles-------------------------------------------
dir_infercnv_output <- "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/snRNA_Processed_Data/InferCNV/outputs/"

genes_loss <- c("VHL", "PBRM1", "BAP1", "SETD2",
                "HIF1A",
                "CDKN2A", 
                "PTEN", 
                "NEGR1",
                "QKI",
                "CADM2", 
                "PTPRD", 
                "NRXN3")
genes_gain <- c("PRKCI", "MECOM",
                "MDM4",
                "MYC",
                "JAK2",
                "SQSTM1",
                "FGFR4")
chr_region <- c(rep("3p", 4),
                "14q",
                "9p21", 
                "10q23",
                "1p31",
                "6q24",
                "3p12",
                "9p23",
                "14q24",
                "3p26", "3p26",
                "1q32",
                "8q24",
                "9q24",
                "5q", "5q")
ccrcc_cna_genes_df <- data.frame(gene_symbol = c(genes_loss, genes_gain),
                                 gene_cna_type = c(rep("Loss", length(genes_loss)), rep("Gain", length(genes_gain))),
                                 chr_region = chr_region)
ccrcc_cna_genes_df

map_bicseq2_log2_copy_ratio2category <- function(log2cr) {
  cnv_cat <- vector(mode = "character", length = length(log2cr))
  cnv_cat[is.na(log2cr)] <- "Not Available"
  cnv_cat[log2cr > -0.05 & log2cr < 0.05] <- "Neutral"
  cnv_cat[log2cr <= -0.6] <- "Deep Loss"
  cnv_cat[log2cr <= -0.05 & log2cr > -0.6] <- "Shallow Loss"
  cnv_cat[log2cr >= 0.05 & log2cr < 0.4] <- "Low Gain"
  cnv_cat[log2cr >= 0.4] <- "High Gain"
  return(cnv_cat)
}
PuBu_colors <- RColorBrewer::brewer.pal(n = 9, name = "PuBu")
PuRd_colors <- RColorBrewer::brewer.pal(n = 9, name = "PuRd")

cna_state_colors <- c("Deep Loss" = PuBu_colors[9],
                      "Shallow Loss" = PuBu_colors[5],
                      "Neutral" = PuBu_colors[3],
                      "Low Gain" = PuRd_colors[5],
                      "High Gain" = PuRd_colors[9],
                      "Not Available" = "grey50")

map_infercnv_state2category <- function(copy_state) {
  cnv_cat <- vector(mode = "character", length = length(copy_state))
  cnv_cat[is.na(copy_state)] <- "Not Available"
  cnv_cat[copy_state == 1] <- "2 Copies"
  cnv_cat[copy_state == 0] <- "0 Copies"
  cnv_cat[copy_state == 0.5] <- "1 Copy"
  cnv_cat[copy_state == 1.5] <- "3 Copies"
  cnv_cat[copy_state == 2] <- "4 Copies"
  cnv_cat[copy_state == 3] <- ">4 Copies"
  return(cnv_cat)
}


copy_number_colors <-  c("0 Copies" = PuBu_colors[9],
                         "1 Copy" = PuBu_colors[5],
                         "2 Copies" = PuBu_colors[3],
                         "3 Copies" = PuRd_colors[5], 
                         "4 Copies" = PuRd_colors[7],
                         ">4 Copies" = PuRd_colors[9],
                         "Not Available" = "grey50")

# map_infercnv_state2category <- function(copy_state) {
#   cnv_cat <- vector(mode = "character", length = length(copy_state))
#   cnv_cat[is.na(copy_state)] <- "Not Available"
#   cnv_cat[copy_state == 1] <- "Neutral"
#   cnv_cat[copy_state == 0] <- "Complete Loss"
#   cnv_cat[copy_state == 0.5] <- "Loss of one copy"
#   cnv_cat[copy_state == 1.5] <- "Addition of one copy"
#   cnv_cat[copy_state == 2] <- "Addition of two copies"
#   cnv_cat[copy_state == 3] <- "Addition > two copies"
#   return(cnv_cat)
# }
# 
# 
# copy_number_colors <-  c("Complete Loss" = PuBu_colors[9],
#                    "Loss of one copy" = PuBu_colors[5],
#                    "Neutral" = PuBu_colors[3],
#                    "Addition of one copy" = PuRd_colors[5], 
#                    "Addition of two copies" = PuRd_colors[7],
#                    "Addition > two copies" = PuRd_colors[9],
#                    "Not Available" = "grey50")




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

get_somatic_mutation_detailed_matrix <- function(pair_tab, maf) {
  genes4mat <- unique(unlist(pair_tab))
  length(genes4mat)
  
  maf <- maf[maf$Hugo_Symbol %in% genes4mat,]
  nrow(maf)
  maf$sampID <- str_split_fixed(string = maf$Tumor_Sample_Barcode, pattern = "_", 2)[,1]
  
  mut_mat <- reshape2::dcast(data = maf, Hugo_Symbol ~ sampID, fun =  function(x) {
    variant_class <- paste0(unique(x), collapse = ",")
    return(variant_class)
  }, value.var = "HGVSp_Short", drop=FALSE)
  rownames(mut_mat) <- as.vector(mut_mat$Hugo_Symbol)
  return(mut_mat)
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


# calculation -------------------------------------------------------------
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


FDR_by_id_columns <- function(p_vector, id_columns, df) {
  ## make a replicate the id columns and make it charactors
  df_id <- matrix(data = as.character(unlist(df[, id_columns])), nrow = nrow(df), ncol = length(id_columns), byrow = F)
  df_id <- data.frame(df_id); colnames(df_id) <- id_columns
  df_id_combo <- data.frame(table(df_id))
  df_id_combo <- df_id_combo[df_id_combo$Freq > 0,]
  ## give a number for each combo
  df_id_combo$combo_id <- 1:nrow(df_id_combo)
  df_id_combo
  df_id <- merge(df_id, df_id_combo[, c(id_columns, "combo_id")], all.x = T)
  if (any(is.na(df_id$combo_id))) {
    stop()
  }
  
  ## for every combo of values in id columns, adjust a set of pvalues
  fdr_vector <- vector(mode = "numeric", length = length(p_vector)) + NA
  for (i in 1:length(unique(df_id$combo_id))) {
    row2adjust <- (df_id$combo_id == i & !is.na(p_vector))
    fdr_vector[row2adjust] <- p.adjust(p = p_vector[row2adjust], method = "fdr")
  }
  return(fdr_vector)
}


