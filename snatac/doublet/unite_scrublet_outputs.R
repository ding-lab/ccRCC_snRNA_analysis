# Yige Wu @WashU March 2020

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
source("./ccRCC_snRNA_analysis/plotting.R")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)


# input dependeciens ------------------------------------------------------
## input scrublet cutoffs used
scrublet_cutoffs_df <- readxl::read_excel(path = "../ccRCC_snATAC/Resources/snATAC_Processed_Data/Doublet_Removal/notes/scrublet_run20201201_cutoffs.xlsx")
## specify input directories
dir_scrublet_out <- "../ccRCC_snATAC/Resources/snATAC_Processed_Data/Doublet_Removal/outputs/"


# preprocess --------------------------------------------------------------
scrublet_cutoffs_used_df <- scrublet_cutoffs_df %>%
  filter(used_for_downstream == "Yes" & call_doublet == "Yes")

# input -------------------------------------------------------------------
scrublet_combined_df <- NULL
for (aliquot_tmp in scrublet_cutoffs_used_df$aliquot) {
  cutoff_used <- scrublet_cutoffs_used_df$scrublet_cutoff[scrublet_cutoffs_used_df$aliquot == aliquot_tmp]
  path_scrublet_out <- paste0(dir_scrublet_out, aliquot_tmp, "/cutoff", cutoff_used, "/", aliquot_tmp, "_scrublet_output_table.csv.gz")
  scrublet_df_tmp <- fread(data.table = F, input = path_scrublet_out)
  scrublet_df_tmp$aliquot <- aliquot_tmp
  scrublet_combined_df <- rbind(scrublet_combined_df, scrublet_df_tmp)
}
table(scrublet_combined_df$aliquot)

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "snATAC_Scrublet_merged_output.", run_id, ".tsv")
write.table(x = scrublet_combined_df, file = file2write, quote = F, sep = "\t", row.names = F)
