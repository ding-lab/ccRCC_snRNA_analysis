# Yige Wu @WashU Jul 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# write paths -------------------------------------------------------------
ids_case <- c("C3N-01200", "C3L-00088", "C3N-01213", "", "")
ids_run <- c("C3N-01200", "C3L-00088", "C3N-01213", "3Sample_PT_and_Tumorcells_ByCellTypeSampleID", "3Sample_AllEpiCells_ByCellTypeSampleID")
paths_monocle_objs <- data.frame(Case = ids_case,
                                 Id_Run = ids_run,
                                 Path_Box = paste0(dir_base, 
                                                   "Resources/snRNA_Processed_Data/Monocle/outputs/",
                                                   ids_run,
                                                   "/combined_subset_pseudotime_qval_1e-10.rds"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "Paths_to_Monocle_Objects.", run_id, ".tsv")
write.table(x = paths_monocle_objs, file = file2write, quote = F, sep = "\t", row.names = F)

