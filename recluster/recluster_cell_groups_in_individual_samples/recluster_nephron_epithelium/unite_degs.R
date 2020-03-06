# Yige Wu @WashU March 2020
## fmake the matrix with the log2 fold change 

# set up libraries and output directory -----------------------------------
## set working directory
baseD = "~/Box/"
setwd(baseD)
source("./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/ccRCC_snRNA_analysis/ccRCC_snRNA_shared.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input the paths for the deg tables
dir_deg_files <- "./Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/Resources/Analysis_Results/recluster/recluster_cell_groups_in_individual_samples/recluster_nephron_epithelium/findallmarkers_malignant_nephron_epithelium_cells_reclustered/20200225.v1/"
### select the files to loop
deg_files2process <- list.files(path = dir_deg_files)
deg_files2process
deg_files2process <- deg_files2process[grepl(pattern = "Wilcox", x = deg_files2process)]
deg_files2process


# loop by aliquot ---------------------------------------------------------
## initiate super table
deg_sup_df <- NULL
for (deg_filename in deg_files2process) {
  aliquot_tmp <- str_split(string = deg_filename, pattern = "\\.")[[1]][1]
  aliquot_tmp
  
  ## input deg
  path_deg_file <- paste0(dir_deg_files, deg_filename)
  deg_df <- fread(path_deg_file)
  
  ## filter by adjusted p value
  deg_df <- deg_df %>%
    filter(p_val_adj < 0.05) %>%
    mutate(aliquot = aliquot_tmp)
  
  deg_sup_df <- rbind(deg_sup_df, deg_df)
}

# write table -------------------------------------------------------------
write.table(x = deg_sup_df, file = paste0(dir_out, "FindAllMarkers.Wilcox.Pos.", run_id, ".tsv"), row.names = F, quote = F, sep = "\t")

