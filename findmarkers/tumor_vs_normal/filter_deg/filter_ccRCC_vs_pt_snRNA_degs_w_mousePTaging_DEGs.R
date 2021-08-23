# Yige Wu @WashU Aug 2021

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
deg_ccRCC_vs_pt_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/annotate_deg/annotate_ccRCC_vs_pt_snRNA_degs_bysex/20210820.v1/ccRCC_vs_PT.snRNA_DEGs.BySex.Summary.20210820.v1.tsv")
deg_pt_age_df <- fread(data.table = F, input = "./Resources/Analysis_Results/findmarkers/tumor_vs_normal/annotate_deg/annotate_mouse_PT_aging_snRNA_degs_bysex/20210820.v1/MousePTAging.snRNA_DEGs.BySex.Summary.20210820.v1.tsv")

# preprocess --------------------------------------------------------------
deg_pt_age_df$genesymbol.human <- toupper(x = deg_pt_age_df$genesymbol.mouse)

# get DEGs common to both sexes -------------------------------------------
deg_common_df <- merge(x = deg_ccRCC_vs_pt_df %>%
                         filter(deg_cat_bysex %in% c("Both_Up", "Both_Down")),
                       y = deg_pt_age_df %>%
                         filter(deg_cat_bysex %in% c("Both_Up", "Both_Down")), by = c("genesymbol.human"), suffixes = c(".ccRCCvsPT", ".mousePTaging"))

# get male specific -------------------------------------------------------
deg_malespec_df <- merge(x = deg_ccRCC_vs_pt_df %>%
                         filter(grepl(x = deg_cat_bysex, pattern = "Male")),
                       y = deg_pt_age_df %>%
                         filter(grepl(x = deg_cat_bysex, pattern = "Male")), by = c("genesymbol.human"), suffixes = c(".ccRCCvsPT", ".mousePTaging"))

# get female specific -------------------------------------------------------
deg_femalespec_df <- merge(x = deg_ccRCC_vs_pt_df %>%
                           filter(grepl(x = deg_cat_bysex, pattern = "Female")),
                         y = deg_pt_age_df %>%
                           filter(grepl(x = deg_cat_bysex, pattern = "Female")), by = c("genesymbol.human"), suffixes = c(".ccRCCvsPT", ".mousePTaging"))

# write output ------------------------------------------------------------
file2write <- paste0(dir_out, "SexesCommon.ccRCC_vs_PT.snRNA_DEGs.OverlapPTAgingDEGs.", run_id, ".tsv")
write.table(x = deg_common_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "MaleSpecific.ccRCC_vs_PT.snRNA_DEGs.OverlapPTAgingDEGs.", run_id, ".tsv")
write.table(x = deg_malespec_df, file = file2write, quote = F, sep = "\t", row.names = F)
file2write <- paste0(dir_out, "FemaleSpecific.ccRCC_vs_PT.snRNA_DEGs.OverlapPTAgingDEGs.", run_id, ".tsv")
write.table(x = deg_femalespec_df, file = file2write, quote = F, sep = "\t", row.names = F)