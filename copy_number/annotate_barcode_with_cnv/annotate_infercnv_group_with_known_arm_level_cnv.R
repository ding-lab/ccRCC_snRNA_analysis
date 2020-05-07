# Yige Wu @WashU May 2020
## make long data frame for each infercnv group whether it has gain or loss of a particular chromosome arm

# set up libraries and output directory -----------------------------------
## set working directory
dir_base = "~/Box/Ding_Lab/Projects_Current/RCC/ccRCC_snRNA/"
setwd(dir_base)
source("./ccRCC_snRNA_analysis/load_pkgs.R")
source("./ccRCC_snRNA_analysis/functions.R")
source("./ccRCC_snRNA_analysis/variables.R")
library("D3GB")
## set run id
version_tmp <- 1
run_id <- paste0(format(Sys.Date(), "%Y%m%d") , ".v", version_tmp)
## set output directory
dir_out <- paste0(makeOutDir(), run_id, "/")
dir.create(dir_out)

# input dependencies ------------------------------------------------------
## input transformed CNV cytoband by infercnv group
group2cnv_region_df <- fread(data.table = F, input = "./Resources/Analysis_Results/copy_number/annotate_barcode_with_cnv/transform_cnv_region2cytoband_by_infercnv_group/20200505.v1/CNV_Region2Cytoband_by_InferCNV_Group20200505.v1.tsv")
## input known arm-level CNvs
knowncnvregion_df <- readxl::read_xlsx(path = "./Resources/Known_Genetic_Alterations/Known_CNV.20200505.v1.xlsx", sheet = "Region")

# summarize ---------------------------------------------------------------
armcnvsummary_df <- group2cnv_region_df %>%
  mutate(cnv_state_binary = ifelse(state > 1, "Gain", "Loss")) %>%
  mutate(start_arm = ifelse(grepl(pattern = "p", x = start_cytoband), paste0(gsub(x = chr, pattern = "chr", replacement = ""), "p"), paste0(gsub(x = chr, pattern = "chr", replacement = ""), "q"))) %>%
  mutate(end_arm = ifelse(grepl(pattern = "p", x = end_cytoband), paste0(gsub(x = chr, pattern = "chr", replacement = ""), "p"), paste0(gsub(x = chr, pattern = "chr", replacement = ""), "q"))) %>%
  select(cell_group_name, cnv_state_binary, start_arm, end_arm)
armcnvsummary_df <- melt(armcnvsummary_df, id.vars = c("cell_group_name", "cnv_state_binary"))
armcnvsummary_df <- unique(armcnvsummary_df %>%
                             select(cell_group_name, cnv_state_binary, value))
test <- table(armcnvsummary_df[, c("cnv_state_binary", "value")])
test <- as.data.frame(test)
