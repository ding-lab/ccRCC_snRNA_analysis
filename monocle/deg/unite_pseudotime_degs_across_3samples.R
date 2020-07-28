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

# input dependencies ------------------------------------------------------
ids_run <- c("C3N-01200", "C3L-00088", "C3N-01213")

# input deg file by case and overlap --------------------------------------
deg_united_df <- NULL
for (id_run in ids_run) {
  path_deg <- paste0("./Resources/snRNA_Processed_Data/Monocle/outputs/", id_run, "/pseudotime_associated_genes/DE_genes_ModelBy_Pseudotime.txt")
  deg_df <- fread(data.table = F, path_deg)
  deg_df <- deg_df %>%
    select(-V1) %>%
    mutate(Id_Run = id_run)
  deg_united_df <- rbind(deg_united_df, deg_df)
}

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "PT_and_TumorCells_ByCellType.Pseudotim_DE_genes.Combine3Samples.", run_id, ".tsv")
write.table(x = deg_united_df, file = file2write, sep = "\t", quote = F, row.names = F)

