# Yige Wu @WashU Sep 2020

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

# input dependencies ------------------------------------------------------
## input the complex annotation from omnipath
complex_omnipath_df <- fread(data.table = F, input = "./Resources/Knowledge/PPI/complexes.txt")

# get the names in cellphone ----------------------------------------------
complex_cp_df <- complex_omnipath_df %>%
  filter(grepl(x = identifiers, pattern = "CellPhoneDB"))
complex_cp_df$identifiers_CellPhoneDB <- sapply(X = complex_cp_df$identifiers, FUN = function(x){
  id_texts_vec <- str_split(string = x, pattern = ";")[[1]]
  idx_cp <- grepl(pattern = "CellPhone", x = id_texts_vec)
  id_text_cp <- id_texts_vec[idx_cp]
  id_text_cp_sim <- gsub(x = id_text_cp, pattern = "CellPhoneDB:", replacement = "")
  return(id_text_cp_sim)
})

# make a data frame with one gene symbol per row --------------------------


# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "CellPhoneDB_complexs_genesymbols.tsv")
write.table(x = complex_cp_df, file = file2write, sep = "\t", quote = F, row.names = F)

