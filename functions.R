# Yige Wu @ WashU 2020 Apr

# make output directory ---------------------------------------------------
makeOutDir_katmai = function(path_script) {
  folders <- strsplit(x = path_script, split = "\\/")[[1]]
  folder_num <- which(folders == "ccRCC_snRNA_analysis") + 1
  dir_analysis_resultnow <- paste(strsplit(paste(folders[folder_num:length(folders)], collapse = "/"), split = "\\.")[[1]][1], sep = "/")
  dir_analysis_resultnow <- paste0(dir_analysis_result, dir_analysis_resultnow, "/")
  dir.create(dir_analysis_resultnow)
  dir_analysis_resultnow_son <- dir_analysis_resultnow
  dirs2make <- NULL
  while (!dir.exists(dir_analysis_resultnow_son)) {
    tmp <- strsplit(dir_analysis_resultnow_son, split = "\\/")[[1]]
    dir_analysis_resultnow_parent <-paste(tmp[-length(tmp)], collapse = "/")
    dir.create(dir_analysis_resultnow_parent)
    dir.create(dir_analysis_resultnow_son)
    dir.create(dir_analysis_resultnow)
    if (!dir.exists(dir_analysis_resultnow_son)) {
      dirs2make[length(dirs2make) + 1] <- dir_analysis_resultnow_son
    }
    dir_analysis_resultnow_son <- dir_analysis_resultnow_parent
  }
  
  if (length(dirs2make) > 0){
    for (i in 1:length(dirs2make)) {
      dir.create(dirs2make[i])
    }
  } 
  return(dir_analysis_resultnow)
}

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
