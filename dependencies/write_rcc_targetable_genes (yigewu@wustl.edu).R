# Yige Wu @WashU Aug 2020
## write druggable targets
## source: https://www.nature.com/articles/nrneph.2017.82#f1

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

# make table --------------------------------------------------------------
targeted_therapy_df <- data.frame(gene_symbol = c("HGF", "MET",
                                                  "GAG6", "AXL",
                                                  "GDNF", "RET",
                                                  "VEGFA", "VEGFB", "VEGFC", "FLT1", "FLT2", "FLT3"), 
                                  ligandreceptor_type = c("ligand", "receptor"),
                                  ligandreceptor_partner = c("MET", "HGF"),
                                  therapy_category = "Targeted_Therapy")


# write output -------------------------------------------------------------
write.table(x = hif_targets, file = paste0(dir_out, "HIF_Target_Genes.", run_id, ".tsv"), quote = F, sep = "\t", row.names = F)


