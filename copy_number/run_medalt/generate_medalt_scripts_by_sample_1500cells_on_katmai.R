# Yige Wu @WashU Jun 2020

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

# input barcode2celltype --------------------------------------------------
barcode2celltype_df <- fread(data.table = F, input = "./Resources/Analysis_Results/annotate_barcode/map_celltype_to_all_cells/20200626.v1/30AliquotIntegration.Barcode2CellType.TumorManualCluster.20200626.v1.tsv")
## get all the downloaded desired output
dir_files <- "./Resources/snRNA_Processed_Data/InferCNV/outputs/Individual.20200305.v1/"
names_files <- list.files(path = dir_files, pattern = "expr.infercnv.dat", recursive = T)
## input id meta data
idmetadata_df <- fread(data.table = F, input = "./Resources/Analysis_Results/sample_info/make_meta_data/20200505.v1/meta_data.20200505.v1.tsv")
## set output directory for medalt
dir_medalt_out <- paste0(dir_base, "Resources/snRNA_Processed_Data/MEDALT/outputs/Individual.20200305.v1.1500tumorcells/")
dir.create(dir_medalt_out)

# extract the tumor only CNVs for each sample -----------------------------
for (name_file in names_files) {
  ## get aliquot id
  id_aliquot <- gsub(x = name_file, pattern = "/expr.infercnv.dat", replacement = "")
  id_aliquot_wu <- idmetadata_df$Aliquot.snRNA.WU[idmetadata_df$Aliquot.snRNA == id_aliquot]
  
  ## make medalt output directory
  dir_medalt_out_aliquot <- paste0(dir_medalt_out, id_aliquot_wu, "/")
  dir.create(dir_medalt_out_aliquot)
  
  ## write script
  file2write <- paste0(dir_out, id_aliquot_wu, ".medalt.sh")
  sink(file = file2write)
  cat(paste0("docker run -v /diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/:/diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/ yigewu/medalt:1.1 ",
             "python /diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/snRNA_Processed_Data/Tree_Analysis/scripts/MEDALT/scTree_setrecursionlimit.py ",
             "-P /diskmnt/Projects/ccRCC_scratch/ccRCC_snRNA/Resources/snRNA_Processed_Data/Tree_Analysis/scripts/MEDALT/ ",
             "-I ", dir_base, "Resources/Analysis_Results/copy_number/run_medalt/generate_medalt_input_files_by_sample_1500cells_on_katmai/", "20200708.v1/",  id_aliquot_wu, ".1500.tumorcells.expr.infercnv.dat ",
             "-O ", dir_medalt_out_aliquot, " -D R -G hg38 ",
             "&> ", dir_medalt_out_aliquot, "Log.$(date +%Y%m%d%H%M%S).txt&"))
  sink()
}

